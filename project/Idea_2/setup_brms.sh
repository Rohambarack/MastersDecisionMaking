#setup c++ toolchain for Rstan on ucloud
sudo add-apt-repository ppa:marutter/rrutter4.0
sudo add-apt-repository ppa:c2d4u.team/c2d4u4.0+
  sudo apt update
sudo apt install r-cran-rstan
#install required R packages
Rscript setupbrms.R