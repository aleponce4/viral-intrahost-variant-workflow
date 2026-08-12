# Third-Party Notices

This repository's own source code is licensed under the MIT License (see [`LICENSE`](LICENSE)).

It additionally **vendors** third-party code that is licensed separately. Vendored code is
**not** covered by this project's MIT License and is redistributed under its own terms. Those
terms continue to apply to the vendored files, to copies of them, and to derivative works of
them, regardless of this project's license.

---

## SNPGenie (`assets/snpgenie/`)

| | |
|---|---|
| **Component** | SNPGenie (`snpgenie.pl`) |
| **Author** | Chase W. Nelson |
| **Copyright** | Copyright (C) 2015, 2016, 2017, 2018, 2019 Chase W. Nelson |
| **License** | GNU General Public License, version 3 or later (GPL-3.0-or-later) |
| **Upstream** | https://github.com/chasewnelson/snpgenie |
| **Pinned commit** | `71584c6c9a30b2c159844f210d70eb89df0f4e19` (branch `master`) |
| **Vendored path** | `assets/snpgenie/snpgenie.pl` |
| **Integrity** | SHA-256 recorded in `assets/snpgenie/checksums.sha256` |

`assets/snpgenie/snpgenie.pl` is an unmodified copy of upstream SNPGenie at the commit pinned
in `assets/snpgenie/PINNED_COMMIT`. It is retained here under the **GNU General Public License
version 3** — its own upstream license — and not under this project's MIT License. The full
license text is available at <https://www.gnu.org/licenses/gpl-3.0.txt>, and the license header
and warranty disclaimer remain in place at the top of `snpgenie.pl`.

The workflow invokes `snpgenie.pl` as an external program (`modules/local/snpgenie/run`); it is
executed, not linked into or copied from by this project's own code. This project's own Perl,
Python, R, and Nextflow sources remain MIT-licensed. If you redistribute this repository, or any
subset of it that includes `assets/snpgenie/`, you must comply with the GPLv3 terms for that
directory — including passing on the license text and the corresponding source.

### Citation

If you use the selection-analysis (`SELECTION`) subworkflow, cite SNPGenie, not just this
pipeline:

> Nelson CW, Moncla LH, Hughes AL (2015). SNPGenie: estimating evolutionary parameters to
> detect natural selection using pooled next-generation sequencing data.
> *Bioinformatics* **31**(22):3709–3711. doi:[10.1093/bioinformatics/btv449](https://doi.org/10.1093/bioinformatics/btv449)

BibTeX:

```bibtex
@article{nelson2015snpgenie,
  author  = {Nelson, Chase W. and Moncla, Louise H. and Hughes, Austin L.},
  title   = {SNPGenie: estimating evolutionary parameters to detect natural selection
             using pooled next-generation sequencing data},
  journal = {Bioinformatics},
  year    = {2015},
  volume  = {31},
  number  = {22},
  pages   = {3709--3711},
  doi     = {10.1093/bioinformatics/btv449}
}
```

---

## Containerized tools

All other analysis tools (LoFreq, iVar, samtools, bcftools, BWA, fastp, FastQC, CliqueSNV,
VILOCA, MultiQC, Quarto, and the Python/R scientific stack) are **not** vendored into this
repository. They are pulled at runtime as pinned Biocontainers images declared in
`conf/containers.config`, and each remains under its own upstream license. Nothing in this
repository redistributes their source or binaries.

To list the exact tool versions used by a given run, see `pipeline_info/software_versions.yml`
in that run's output directory, or the Software Provenance section of the generated executive
HTML report.
