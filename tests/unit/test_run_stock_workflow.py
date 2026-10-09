"""Unit tests for scripts/run_stock_workflow.sh.

Run from the repository root:  python -m unittest discover -s tests/unit -v

The launcher chains three Nextflow runs. These tests replace Nextflow with a fake
that records every call and writes the files a real run would, so the order, the
arguments and the hand-off between steps are checked without Nextflow or data.
"""
import gzip
import os
import shutil
import stat
import subprocess
import tempfile
import unittest

REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..'))
SCRIPT = os.path.join(REPO, 'scripts', 'run_stock_workflow.sh')
STAGE_A_PARAMS = os.path.join(REPO, 'assets', 'stage_a_viralmetagenome.params.yaml')

FAKE_NEXTFLOW = r"""#!/usr/bin/env bash
tag="@TAG@"
if [ "$1" = "-version" ]; then echo "      version 99.${tag}.0 build 1 (java ${NXF_JAVA_HOME:-default})"; exit 0; fi
printf '%s|%s|%s|%s\n' "$tag" "$PWD" "${NXF_JAVA_HOME:-unset}" "$*" >> "$FAKE_LOG"
args=("$@")
outdir=""; stock=""; primer=""; input=""
for ((i = 0; i < ${#args[@]}; i++)); do
    case "${args[$i]}" in
        --outdir) outdir="${args[$((i + 1))]}" ;;
        --stock_name) stock="${args[$((i + 1))]}" ;;
        --stock_primer_bed) primer=1 ;;
        --input) input="${args[$((i + 1))]}" ;;
    esac
done
case "$*" in
    *nf-core/viralmetagenome*)
        sample=$(sed -n 2p "$input" | cut -d, -f1)
        mkdir -p "$outdir/consensus/seq/variant-calling/$sample"
        for n in ${FAKE_CLUSTERS-1}; do
            : > "$outdir/consensus/seq/variant-calling/$sample/${sample}_cl${n}_itvariant-calling.consensus.fasta"
        done ;;
    *BUILD_STOCK_REFERENCE*)
        d="$outdir/StockReference/$stock"; mkdir -p "$d/qc"
        : > "$d/$stock.fasta"; : > "$d/$stock.gff3"; : > "$d/qc/$stock.liftover.tsv"
        [ -n "$primer" ] && : > "$d/$stock.primers.bed"
        true ;;
    *) mkdir -p "$outdir" ;;
esac
"""


@unittest.skipUnless(shutil.which('bash'), 'the launcher is a bash script')
class LauncherTests(unittest.TestCase):
    def setUp(self):
        self._tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self._tmp.cleanup)
        self.dir = os.path.realpath(self._tmp.name)
        self.fake_a = self.make_fake('A')
        self.fake_b = self.make_fake('B')
        for name in ('R1.fq.gz', 'R2.fq.gz', 'pool.fa', 'lab.fa', 'lab.gff3',
                     'samples.csv', 'primers.bed', 'my.consensus.fasta', 'extra.config'):
            open(self.p(name), 'w').close()
        self.log = self.p('calls.log')
        self.out = self.p('out')

    def p(self, *parts):
        return os.path.join(self.dir, *parts)

    def make_fake(self, tag):
        path = self.p(f'fake_nextflow_{tag}')
        with open(path, 'w', encoding='utf-8', newline='\n') as f:
            f.write(FAKE_NEXTFLOW.replace('@TAG@', tag))
        os.chmod(path, os.stat(path).st_mode | stat.S_IXUSR)
        return path

    def run_launcher(self, *args, env=None):
        e = dict(os.environ, NXF_STAGE_A=self.fake_a, NXF_STAGE_B=self.fake_b,
                 FAKE_LOG=self.log)
        e.pop('NXF_JAVA_HOME', None)
        e.update(env or {})
        return subprocess.run(['bash', SCRIPT, *args], capture_output=True, text=True,
                              env=e, cwd=self.dir)

    def full(self, *extra):
        return ['--stock', 's1', '--outdir', self.out,
                '--stock-reads', self.p('R1.fq.gz'), self.p('R2.fq.gz'),
                '--reference-pool', self.p('pool.fa'),
                '--lab-fasta', self.p('lab.fa'), '--lab-gff', self.p('lab.gff3'),
                '--samples', self.p('samples.csv'), *extra]

    def calls(self):
        if not os.path.exists(self.log):
            return []
        with open(self.log, encoding='utf-8') as f:
            rows = [line.rstrip('\n').split('|', 3) for line in f if line.strip()]
        return [{'tag': t, 'cwd': c, 'java': j, 'args': a} for t, c, j, a in rows]

    def stock_ref(self, *parts):
        return os.path.join(self.out, 'reference', 'results', 'StockReference', 's1', *parts)

    # --- the chain ---------------------------------------------------------

    def test_three_steps_run_in_order(self):
        result = self.run_launcher(*self.full())
        self.assertEqual(result.returncode, 0, result.stderr)
        calls = self.calls()
        self.assertEqual([c['tag'] for c in calls], ['A', 'B', 'B'])
        self.assertIn('nf-core/viralmetagenome', calls[0]['args'])
        self.assertIn('-entry BUILD_STOCK_REFERENCE', calls[1]['args'])
        self.assertNotIn('-entry', calls[2]['args'])

    def test_stage_a_arguments(self):
        self.run_launcher(*self.full())
        a = self.calls()[0]['args']
        self.assertIn('-r 1.2.0', a)
        self.assertIn('-profile docker', a)
        self.assertIn(f'-params-file {STAGE_A_PARAMS}', a)
        self.assertIn(f'--reference_pool {self.p("pool.fa")}', a)
        self.assertIn(f'--input {self.out}/stage_a/samplesheet.csv', a)

    def test_stage_a_samplesheet_names_the_stock(self):
        self.run_launcher(*self.full())
        with open(os.path.join(self.out, 'stage_a', 'samplesheet.csv'), encoding='utf-8') as f:
            self.assertEqual(f.read(),
                             f'sample,fastq_1,fastq_2\ns1,{self.p("R1.fq.gz")},{self.p("R2.fq.gz")}\n')

    def test_consensus_found_by_stage_a_is_handed_to_the_reference_step(self):
        self.run_launcher(*self.full())
        expected = (f'{self.out}/stage_a/results/consensus/seq/variant-calling/s1/'
                    's1_cl1_itvariant-calling.consensus.fasta')
        ref = self.calls()[1]['args']
        self.assertIn(f'--stock_consensus {expected}', ref)
        self.assertIn('--stock_name s1', ref)
        self.assertIn(f'--fasta {self.p("lab.fa")}', ref)
        self.assertIn(f'--gff {self.p("lab.gff3")}', ref)

    def test_reference_outputs_are_handed_to_the_variant_run(self):
        self.run_launcher(*self.full())
        var = self.calls()[2]['args']
        self.assertIn(f'--fasta {self.stock_ref("s1.fasta")}', var)
        self.assertIn(f'--gff {self.stock_ref("s1.gff3")}', var)
        self.assertIn('--viral_contig s1', var)
        self.assertIn(f'--liftover_tsv {self.stock_ref("qc", "s1.liftover.tsv")}', var)
        self.assertIn(f'--input {self.p("samples.csv")}', var)
        self.assertNotIn('--protocol', var)

    def test_each_step_runs_in_its_own_directory(self):
        self.run_launcher(*self.full())
        self.assertEqual([os.path.basename(c['cwd']) for c in self.calls()],
                         ['stage_a', 'reference', 'variants'])

    def test_every_run_resumes(self):
        self.run_launcher(*self.full())
        self.assertTrue(all('-resume' in c['args'] for c in self.calls()))

    # --- options -----------------------------------------------------------

    def test_primer_bed_moves_onto_the_stock_and_switches_to_amplicon(self):
        self.run_launcher(*self.full('--primer-bed', self.p('primers.bed')))
        calls = self.calls()
        self.assertIn(f'--stock_primer_bed {self.p("primers.bed")}', calls[1]['args'])
        self.assertIn(f'--protocol amplicon --primer_bed {self.stock_ref("s1.primers.bed")}',
                      calls[2]['args'])

    def test_supplied_consensus_skips_stage_a(self):
        args = ['--stock', 's1', '--outdir', self.out, '--consensus', self.p('my.consensus.fasta'),
                '--lab-fasta', self.p('lab.fa'), '--lab-gff', self.p('lab.gff3'),
                '--samples', self.p('samples.csv')]
        result = self.run_launcher(*args)
        self.assertEqual(result.returncode, 0, result.stderr)
        calls = self.calls()
        self.assertEqual([c['tag'] for c in calls], ['B', 'B'])
        self.assertIn(f'--stock_consensus {self.p("my.consensus.fasta")}', calls[0]['args'])
        self.assertNotIn('viralmetagenome', ' '.join(c['args'] for c in calls))

    def test_the_two_stages_can_use_different_nextflow_binaries(self):
        self.run_launcher(*self.full())
        self.assertEqual([c['tag'] for c in self.calls()], ['A', 'B', 'B'])

    def test_extra_arguments_reach_only_their_own_step(self):
        self.run_launcher(*self.full('--stage-a-extra', '--max_cpus 4',
                                     '--variants-extra', '--lofreq_sb_thresh 0'))
        a, ref, var = (c['args'] for c in self.calls())
        self.assertIn('--max_cpus 4', a)
        self.assertNotIn('--lofreq_sb_thresh', a + ref)
        self.assertIn('--lofreq_sb_thresh 0', var)

    def test_config_and_profile_go_to_every_step(self):
        self.run_launcher(*self.full('--config', self.p('extra.config'), '--profile', 'singularity'))
        for c in self.calls():
            self.assertIn(f'-c {self.p("extra.config")}', c['args'])
            self.assertIn('-profile singularity', c['args'])

    def test_site_default_profile_comes_from_the_environment(self):
        self.run_launcher(*self.full(), env={'STOCK_WORKFLOW_PROFILE': 'apptainer'})
        self.assertTrue(all('-profile apptainer' in c['args'] for c in self.calls()))

    def test_profile_flag_beats_the_environment(self):
        self.run_launcher(*self.full('--profile', 'singularity'),
                          env={'STOCK_WORKFLOW_PROFILE': 'apptainer'})
        self.assertTrue(all('-profile singularity' in c['args'] for c in self.calls()))

    def test_stage_a_can_have_its_own_java(self):
        self.run_launcher(*self.full(), env={'NXF_STAGE_A_JAVA_HOME': '/opt/jdk17'})
        a, ref, var = self.calls()
        self.assertEqual(a['java'], '/opt/jdk17')
        self.assertEqual((ref['java'], var['java']), ('unset', 'unset'))

    def test_the_other_steps_can_have_their_own_java(self):
        self.run_launcher(*self.full(), env={'NXF_STAGE_B_JAVA_HOME': '/opt/jdk21'})
        a, ref, var = self.calls()
        self.assertEqual(a['java'], 'unset')
        self.assertEqual((ref['java'], var['java']), ('/opt/jdk21', '/opt/jdk21'))

    def test_the_log_records_which_java_each_stage_used(self):
        self.run_launcher(*self.full(), env={'NXF_STAGE_A_JAVA_HOME': '/opt/jdk17'})
        with open(os.path.join(self.out, 'stock_run.log'), encoding='utf-8') as f:
            text = f.read()
        self.assertIn('Stage A nextflow: version 99.A.0 build 1 (java /opt/jdk17)', text)
        self.assertIn('Stage B nextflow: version 99.B.0 build 1 (java default)', text)

    def test_a_dry_run_shows_the_java_setting(self):
        result = self.run_launcher(*self.full('--dry-run'), env={'NXF_STAGE_A_JAVA_HOME': '/opt/jdk17'})
        self.assertIn('env NXF_JAVA_HOME=/opt/jdk17 ', result.stdout)

    def test_a_single_step_runs_alone(self):
        result = self.run_launcher(*self.full('--step', 'stage-a'))
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual([c['tag'] for c in self.calls()], ['A'])

    def test_reference_step_alone_finds_an_earlier_stage_a(self):
        self.run_launcher(*self.full('--step', 'stage-a'))
        result = self.run_launcher(*self.full('--step', 'reference'))
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual([c['tag'] for c in self.calls()], ['A', 'B'])

    # --- dry run and the log ----------------------------------------------

    def test_dry_run_prints_the_commands_and_writes_nothing(self):
        result = self.run_launcher(*self.full('--dry-run'))
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(self.calls(), [])
        self.assertFalse(os.path.exists(self.out))
        for step in ('Stage A: de novo consensus', 'Stock reference', 'Variant calling'):
            self.assertIn(f'[dry-run] {step}: cd ', result.stdout)
        self.assertIn('nf-core/viralmetagenome', result.stdout)
        self.assertIn('BUILD_STOCK_REFERENCE', result.stdout)
        self.assertIn('CONSENSUS_FROM_STAGE_A', result.stdout)
        self.assertIn('sample,fastq_1,fastq_2', result.stdout)

    def test_log_records_versions_and_commands(self):
        self.run_launcher(*self.full())
        with open(os.path.join(self.out, 'stock_run.log'), encoding='utf-8') as f:
            text = f.read()
        self.assertIn('run started', text)
        self.assertIn('stock s1, step all', text)
        self.assertIn('pipeline commit:', text)
        self.assertIn('Stage A nextflow: version 99.A.0', text)
        self.assertIn('Stage B nextflow: version 99.B.0', text)
        self.assertIn('nf-core/viralmetagenome', text)
        self.assertIn('BUILD_STOCK_REFERENCE', text)

    # --- failures ----------------------------------------------------------

    def test_no_arguments_prints_usage_and_exits_2(self):
        result = self.run_launcher()
        self.assertEqual(result.returncode, 2)
        self.assertIn('Usage:', result.stderr)

    def test_every_real_run_adds_a_header_to_the_log_even_without_stage_a(self):
        self.run_launcher(*self.full('--step', 'stage-a'))
        self.run_launcher(*self.full('--step', 'reference'))
        with open(os.path.join(self.out, 'stock_run.log'), encoding='utf-8') as f:
            text = f.read()
        self.assertEqual(text.count('run started'), 2)
        self.assertIn('step reference', text)

    def test_a_stock_name_that_cannot_be_a_contig_name_is_refused(self):
        for bad in ('my stock', 'a/b', 'x;y', 'tc83*'):
            args = self.full()
            args[args.index('--stock') + 1] = bad
            result = self.run_launcher(*args)
            self.assertEqual(result.returncode, 2, bad)
            self.assertIn('--stock may hold only', result.stderr)
        self.assertEqual(self.calls(), [])

    def test_missing_input_file_is_named(self):
        args = self.full()
        args[args.index('--reference-pool') + 1] = self.p('nope.fa')
        result = self.run_launcher(*args)
        self.assertEqual(result.returncode, 2)
        self.assertIn('--reference-pool', result.stderr)
        self.assertEqual(self.calls(), [])

    def test_variants_step_without_a_reference_says_so(self):
        result = self.run_launcher('--stock', 's1', '--outdir', self.out,
                                   '--samples', self.p('samples.csv'), '--step', 'variants')
        self.assertEqual(result.returncode, 2)
        self.assertIn('--step reference', result.stderr)
        self.assertEqual(self.calls(), [])

    def test_two_consensus_files_stop_the_run_and_are_listed(self):
        result = self.run_launcher(*self.full(), env={'FAKE_CLUSTERS': '1 2'})
        self.assertEqual(result.returncode, 2)
        self.assertIn('more than one consensus', result.stderr)
        self.assertIn('s1_cl1_itvariant-calling.consensus.fasta', result.stderr)
        self.assertIn('s1_cl2_itvariant-calling.consensus.fasta', result.stderr)
        self.assertEqual([c['tag'] for c in self.calls()], ['A'])

    def test_no_consensus_stops_the_run(self):
        result = self.run_launcher(*self.full(), env={'FAKE_CLUSTERS': ''})
        self.assertEqual(result.returncode, 2)
        self.assertIn('no consensus found', result.stderr)
        self.assertEqual([c['tag'] for c in self.calls()], ['A'])

    def test_consensus_with_the_stage_a_step_is_refused(self):
        result = self.run_launcher(*self.full('--step', 'stage-a', '--consensus', self.p('my.consensus.fasta')))
        self.assertEqual(result.returncode, 2)
        self.assertIn('nothing to do', result.stderr)

    def test_a_primer_bed_without_a_primer_reference_is_caught_before_the_variant_run(self):
        self.run_launcher(*self.full())
        result = self.run_launcher(*self.full('--step', 'variants', '--primer-bed', self.p('primers.bed')))
        self.assertEqual(result.returncode, 2)
        self.assertIn('without --primer-bed', result.stderr)
        self.assertEqual([c['tag'] for c in self.calls()], ['A', 'B', 'B'])

    def test_a_failing_step_stops_the_chain(self):
        failing = self.p('fail_nextflow')
        with open(failing, 'w', encoding='utf-8', newline='\n') as f:
            f.write('#!/usr/bin/env bash\n[ "$1" = "-version" ] && exit 0\nexit 7\n')
        os.chmod(failing, 0o755)
        result = self.run_launcher(*self.full(), env={'NXF_STAGE_A': failing})
        self.assertNotEqual(result.returncode, 0)
        self.assertEqual(self.calls(), [])


    # --- Stage A subsample and per-stage configs ---------------------------

    def write_fastq_gz(self, name, reads, tag):
        """A real gzip FASTQ with `reads` records, so --stage-a-pairs has something to cut."""
        path = self.p(name)
        with gzip.open(path, 'wt', encoding='utf-8') as f:
            for i in range(reads):
                f.write(f'@{tag}_read{i}\nACGTACGT\n+\nIIIIIIII\n')
        return path

    def read_ids(self, path):
        with gzip.open(path, 'rt', encoding='utf-8') as f:
            return [line.strip() for n, line in enumerate(f) if n % 4 == 0]

    def with_real_reads(self, reads=10):
        r1 = self.write_fastq_gz('real_R1.fastq.gz', reads, 'p')
        r2 = self.write_fastq_gz('real_R2.fastq.gz', reads, 'p')
        args = self.full()
        i = args.index('--stock-reads')
        args[i + 1], args[i + 2] = r1, r2
        return args

    def test_stage_a_pairs_keeps_the_first_pairs_and_points_stage_a_at_them(self):
        result = self.run_launcher(*self.with_real_reads(10), '--stage-a-pairs', '4')
        self.assertEqual(result.returncode, 0, result.stderr)
        sub = os.path.join(self.out, 'stage_a', 'subsample')
        for r in ('R1', 'R2'):
            self.assertEqual(self.read_ids(os.path.join(sub, f's1_{r}.fastq.gz')),
                             ['@p_read0', '@p_read1', '@p_read2', '@p_read3'])
        with open(os.path.join(self.out, 'stage_a', 'samplesheet.csv'), encoding='utf-8') as f:
            self.assertIn(f's1,{sub}/s1_R1.fastq.gz,{sub}/s1_R2.fastq.gz', f.read())

    def test_stage_a_pairs_does_not_touch_the_variant_samples(self):
        self.run_launcher(*self.with_real_reads(10), '--stage-a-pairs', '4')
        var = self.calls()[2]['args']
        self.assertIn(f'--input {self.p("samples.csv")}', var)

    def test_stage_a_pairs_larger_than_the_file_uses_all_reads_and_says_so(self):
        result = self.run_launcher(*self.with_real_reads(10), '--stage-a-pairs', '100')
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn('only 10 read pairs', result.stdout)
        sub = os.path.join(self.out, 'stage_a', 'subsample', 's1_R1.fastq.gz')
        self.assertEqual(len(self.read_ids(sub)), 10)

    def test_stage_a_pairs_is_recorded_in_the_log(self):
        self.run_launcher(*self.with_real_reads(10), '--stage-a-pairs', '4')
        with open(os.path.join(self.out, 'stock_run.log'), encoding='utf-8') as f:
            self.assertIn('Stage A reads: the first 4 read pairs of', f.read())

    def test_stage_a_pairs_is_reused_when_nothing_changed(self):
        args = self.with_real_reads(10)
        self.run_launcher(*args, '--stage-a-pairs', '4')
        target = os.path.join(self.out, 'stage_a', 'subsample', 's1_R1.fastq.gz')
        before = os.stat(target).st_mtime_ns
        result = self.run_launcher(*args, '--stage-a-pairs', '4')
        self.assertIn('reusing the first 4 read pairs', result.stdout)
        self.assertEqual(os.stat(target).st_mtime_ns, before)

    def test_stage_a_pairs_is_redone_when_the_number_changes(self):
        args = self.with_real_reads(10)
        self.run_launcher(*args, '--stage-a-pairs', '4')
        self.run_launcher(*args, '--stage-a-pairs', '6')
        target = os.path.join(self.out, 'stage_a', 'subsample', 's1_R1.fastq.gz')
        self.assertEqual(len(self.read_ids(target)), 6)

    def test_stage_a_pairs_values_that_are_not_positive_whole_numbers_are_refused(self):
        for bad in ('0', 'abc', '-3', '2.5', ''):
            result = self.run_launcher(*self.full(), '--stage-a-pairs', bad)
            self.assertEqual(result.returncode, 2, repr(bad))
            self.assertIn('positive whole number', result.stderr)
        self.assertEqual(self.calls(), [])

    def test_stage_a_pairs_with_a_leading_zero_is_read_as_decimal(self):
        result = self.run_launcher(*self.with_real_reads(10), '--stage-a-pairs', '08')
        self.assertEqual(result.returncode, 0, result.stderr)
        sub = os.path.join(self.out, 'stage_a', 'subsample', 's1_R1.fastq.gz')
        self.assertEqual(len(self.read_ids(sub)), 8)

    def test_stage_a_pairs_with_a_file_that_is_not_gzip_is_refused(self):
        result = self.run_launcher(*self.full(), '--stage-a-pairs', '4')
        self.assertEqual(result.returncode, 2)
        self.assertIn('as a gzip FASTQ', result.stderr)
        self.assertEqual(self.calls(), [])

    def test_stage_a_pairs_with_r1_and_r2_of_different_length_is_refused(self):
        args = self.with_real_reads(10)
        i = args.index('--stock-reads')
        args[i + 2] = self.write_fastq_gz('short_R2.fastq.gz', 3, 'p')
        result = self.run_launcher(*args, '--stage-a-pairs', '5')
        self.assertEqual(result.returncode, 2)
        self.assertIn('different numbers of reads', result.stderr)

    def test_stage_a_pairs_is_ignored_when_a_consensus_replaces_stage_a(self):
        args = ['--stock', 's1', '--outdir', self.out, '--consensus', self.p('my.consensus.fasta'),
                '--lab-fasta', self.p('lab.fa'), '--lab-gff', self.p('lab.gff3'),
                '--samples', self.p('samples.csv'), '--stage-a-pairs', '4']
        result = self.run_launcher(*args)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertFalse(os.path.exists(os.path.join(self.out, 'stage_a')))

    def test_dry_run_describes_the_subsample_and_writes_nothing(self):
        result = self.run_launcher(*self.with_real_reads(10), '--stage-a-pairs', '4', '--dry-run')
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn('would keep the first 4 read pairs', result.stdout)
        self.assertIn('s1_R1.fastq.gz', result.stdout)
        self.assertFalse(os.path.exists(self.out))

    def test_per_stage_configs_reach_only_their_own_stage(self):
        for name in ('all.config', 'a.config', 'v.config'):
            open(self.p(name), 'w').close()
        self.run_launcher(*self.full('--config', self.p('all.config'),
                                     '--stage-a-config', self.p('a.config'),
                                     '--variants-config', self.p('v.config')))
        a, ref, var = (c['args'] for c in self.calls())
        self.assertIn(f'-c {self.p("all.config")}', a)
        self.assertIn(f'-c {self.p("a.config")}', a)
        self.assertNotIn('v.config', a)
        self.assertIn(f'-c {self.p("all.config")}', ref)
        self.assertNotIn('a.config', ref)
        self.assertNotIn('v.config', ref)
        self.assertIn(f'-c {self.p("all.config")}', var)
        self.assertIn(f'-c {self.p("v.config")}', var)
        self.assertNotIn('a.config', var)

    def test_the_per_stage_config_comes_after_the_shared_one_so_it_wins(self):
        for name in ('all.config', 'a.config'):
            open(self.p(name), 'w').close()
        self.run_launcher(*self.full('--config', self.p('all.config'),
                                     '--stage-a-config', self.p('a.config')))
        a = self.calls()[0]['args']
        self.assertLess(a.index('all.config'), a.index('a.config'))

    def test_a_missing_per_stage_config_is_named(self):
        for flag in ('--stage-a-config', '--variants-config'):
            result = self.run_launcher(*self.full(flag, self.p('nope.config')))
            self.assertEqual(result.returncode, 2, flag)
            self.assertIn(flag, result.stderr)
        self.assertEqual(self.calls(), [])


if __name__ == '__main__':
    unittest.main()
