python3 -mvenv /tmp/coal_hmm_review
source /tmp/coal_hmm_review/bin/activate
pip3 install -U pip wheel cython
pip3 install -r requirements.txt
pip3 install git+https://github.com/popgenmethods/smcpp@v1.15.3
git clone https://github.com/lh3/psmc /tmp/coal_hmm_review/psmc
make -C /tmp/coal_hmm_review/psmc
# Attempt to acquire MSMC & related tools
mkdir -p /tmp/coal_hmm_review/msmc
wget -o /tmp/coal_hmm_review/msmc/msmc https://github.com/stschiff/msmc/releases/download/v1.1.0/msmc_1.1.0_linux64bit
chmod +x /tmp/coal_hmm_review/msmc/msmc
git clone https://github.com/stschiff/msmc-tools /tmp/coal_hmm_review/msmc/tools
