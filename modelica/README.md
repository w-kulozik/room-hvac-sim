# Modelica cross-validation

`RoomSim.mo` is an independent implementation of the room model (balances, derived parameters
and the sinusoidal summer day) in Modelica, written from [docs/model.md](../docs/model.md).
It is used to cross-validate the Python model (decision D58 in
[docs/decisions.md](../docs/decisions.md); results: D60 and
`figures/modelica_crossvalidation.md`).

- `results/`: the OpenModelica results of the two cases (CSV, every 600 s for three days) and
  the OpenModelica version. They are committed, so the comparison runs without OpenModelica.
- `python scripts/crossvalidate_modelica.py` compares; with `--run-modelica` it first
  regenerates the results (OpenModelica's `omc` from `$OPENMODELICAHOME/bin` or the PATH).
- `tests/test_modelica.py` checks the bounds of D58.
