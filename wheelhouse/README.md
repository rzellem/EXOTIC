# Local UltraNest wheels

These Windows x64 wheels package UltraNest 4.5.0 for the two CPython runtimes
supported by this checkout:

- `ultranest-4.5.0-cp312-cp312-win_amd64.whl` — CPython 3.12
- `ultranest-4.5.0-cp313-cp313-win_amd64.whl` — CPython 3.13

They were built locally from the UltraNest 4.5.0 source distribution with the
Visual Studio 2022 C++ toolchain. The compiled `ultranest.mlfriends` extension
was imported successfully under both matching interpreters.

SHA-256:

- CPython 3.12: `69C57FC84FD475AEA6A959F473F942192A6839860849D626456E0F61406B06FF`
- CPython 3.13: `09CA8F766E039F072B9182B6072F0E318FA6173EA516368B988E5F9F56A9834C`

`run_exotic_windows.bat` adds this directory as a pip `--find-links` source and
requires a binary UltraNest distribution, preventing a fallback source build.
