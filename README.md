# xsat2

## Install Dependencies

```bash
apt install python3-venv -y
python3 -m venv xsat
source xsat/bin/activate
pip install --upgrade pip
pip install -r requirements.txt
```

## Batch Test

```bash
# small benchmarks
python test_benchmarks.py --small
# middle benchmarks
python test_benchmarks.py --middle
# large benchmarks
python test_benchmarks.py --large
```

## Single Benchmark Test

```bash
make IN=Benchmarks/middle/div3.c.50.smt2
python xsat.py 
```
