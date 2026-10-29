# AI-use log

A running record of how Claude Code is used in this project. It is the raw material for the report's AI-use statement, which must say what Claude was asked to do, what was verified independently, and give at least one concrete example where Claude's output had to be corrected or rejected. Add an entry whenever Claude writes a non-trivial piece of code or the owner catches an error.

## Entry template

- **Date:**
- **Asked Claude to:**
- **What Claude produced:**
- **How it was checked** (tests, comparison with the reference function, reading the code against the equations):
- **Problems found and how they were corrected** (leave blank if none):

## Entries

### 2026-10-28: Repository scaffolding

- **Asked Claude to:** read the project description, propose a repository layout, create the folder structure, a `.gitignore`, `CLAUDE.md`, and starter notes, and verify the package name is free on CRAN.
- **What Claude produced:** the folder skeleton, `.gitignore`, `CLAUDE.md`, stub `README.md`, `LICENSE.md`, and the notes in `dev/`. It also ran the name-availability checks (results in `dev/progress_log.md`).
- **How it was checked:** the owner reviewed and approved the plan before any files were created.
- **Problems found and how they were corrected:** the first attempt to install a Python PDF library in a scratch folder was declined by the owner in favor of installing Poppler, which was then used to read the handout. A first installation attempt for `devtools` failed because of locked library files, which led to using a compatible binary version and the per-user library.
