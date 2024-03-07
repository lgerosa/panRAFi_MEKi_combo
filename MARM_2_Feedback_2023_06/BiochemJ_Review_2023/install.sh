conda create -n biochemj
source activate biochemj
conda install python=3.9.12
python setup.py install
conda install ipykernel
python -m ipykernel install --user --name BiochemJ_Review
conda install -y matplotlib cython joblib pandas jinja2 tqdm
conda install -c conda-forge python-libsbml
pip install synergy
conda install -y -c alubbock pysb
conda install -y -c conda-forge ipywidgets
