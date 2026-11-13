# paper_archive

Exact script copies that were used to make the figures and tables of the submitted manuscript
(eLife reviewed preprint 111513). They are kept exactly as they were, including older and
experimental versions, so that every panel can be traced to the code that produced it.

- `scripts/` : one file per distinct script content. Names end in the first 8 characters of the
  file's MD5 checksum (`<name>__<md5>.R`). Two files with the same base name are different versions.
- `script_archive_index.csv` : for each file, its full MD5, how many copies existed in the figure
  folders, when it was last saved, and every folder it was found in.
- `script_archive_record.json` : counts and settings from the run that built this folder.

Scripts were copied byte-for-byte except that your computer's folder names were removed from file paths inside the code
(index columns `had_local_path` and `path_scrubbed`; `md5` is the original file's checksum, `archive_md5` is the archived file's).

The original pipeline is in the numbered folders (00, 01, 02) at the top level of this repository.
These archived copies are for provenance, not for new analyses.
