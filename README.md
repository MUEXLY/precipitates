# precipitates
Repository for computing precipitate sizes in a molecular dynamics run

## Example run on Palmetto 2

Add anaconda module
```bash
module add anaconda3
```

Create virtual environment, activate it, and install dependencies (will take awhile)
```bash
python -m venv venv
source venv/bin/activate
pip install -r requirements.txt
```

Run example (will also take awhile, but has a progress message)
```bash
python clusters.py
```

This example run generates `clusters.pdf`. The run parameters are defined in `clusters.py`.