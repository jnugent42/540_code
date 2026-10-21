#echo "Installing python packages:"
#pip install -r requirements.txt

echo "activating virtual environment:"
source /vols/hyperk/WCTE/WCTE540/lb_git/env/bin/activate

echo "setting up CUDA:"
source /vols/software/cuda/setup.sh 11.2.0
