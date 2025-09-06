# roomsim — dynamic room model for HVAC control

A small, verified simulator of a single room (air node + thermal mass node + moisture balance,
with psychrometrics, outdoor weather, solar gains through a window and occupancy). It is the
plant on which discrete-time air-conditioning controllers (on/off, PI, cascade, MPC) are tested.

Work in progress. See the sections below as they are added.

## Install and test

```bash
pip install -e .[dev]
pytest -v
```
