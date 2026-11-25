# EXOTIC repository instructions

## Keep the Colab wheel current

- Whenever you change EXOTIC Python code, native code, or packaging inputs used
  by Colab, rebuild the Linux x86-64 EXOTIC wheel from that same revision before
  handing off the change. A source-only update leaves Colab on the previous
  wheel. Notebook, documentation, and test-only edits need no rebuild unless
  they also change the installed package.
- Follow the [Colab build instructions](scripts/README.md#prebuilt-colab-wheels).
  Update the EXOTIC wheel, `wheelhouse/colab/index.json`, and
  `wheelhouse/colab/requirements.txt` together. Remove the superseded EXOTIC
  wheel and keep the tested dependency wheels when their versions still satisfy
  the new build.
- Run `python -m pytest -q tests/test_colab_wheelhouse.py` and check that the
  manifest selects the rebuilt wheel and that its installed version identifies
  the source revision. Report the build version and any verification you could
  not run.
- Keep the wheel and manifest updates with the source changes. A rebuild does
  not authorize a push or a release. Push only when the user specifically asks.
