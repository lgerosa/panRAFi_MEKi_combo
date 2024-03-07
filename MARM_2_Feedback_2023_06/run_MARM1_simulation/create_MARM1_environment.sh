module load conda2
conda create -n marm1 -y -c alubbock pysb cython pandas
eval $(conda shell.bash activate marm1)
python -m pip install git+git://github.com/jmuhlich/pysb.git@marm1-integration
