$e="$env:USERPROFILE\miniconda3\envs\plip"
$env:BABEL_LIBDIR="$e\Library\bin"
$env:BABEL_DATADIR="$e\Library\share\openbabel\3.2.1"
$env:PATH="$e;$e\Library\bin;$env:PATH"
Set-Location "$env:TEMP\snaclex-bench\snaclex-main"
& "$e\python.exe" plip_compare.py 2 *> plip_run.log
