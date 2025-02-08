#!/bin/bash -ex
echo "install OS updates, python, nss-tools" >> log.txt
yum update -y
yum install nss-tools -y
yum install python3.11.x86_64 -y
yum install nginx -y
yum install systemctl -y
echo "create venv" >> log.txt
python3.11 -m venv venv
source venv/bin/activate
python3.11 -m ensurepip
echo "download app from github" >> log.txt
curl -O https://raw.githubusercontent.com/bmorris3/jdaviz-aws-demo/refs/heads/main/requirements.txt
curl -O https://raw.githubusercontent.com/bmorris3/jdaviz-aws-demo/refs/heads/main/solara.css
curl -O https://raw.githubusercontent.com/bmorris3/jdaviz-aws-demo/refs/heads/main/app.py
curl -O https://github.com/bmorris3/jdaviz-aws-demo/raw/refs/heads/main/nginx.conf
curl -O https://github.com/bmorris3/jdaviz-aws-demo/raw/refs/heads/main/launch-server.sh
mv nginx.conf /etc/nginx/.
systemctl start nginx.service
python3.11 -m pip install -U pip
echo "install python deps" >> log.txt
python3.11 -m pip install -r requirements.txt
echo "launch server" >> log.txt
./launch-server.sh
