#!/bin/sh
set -ex
wget https://github.com/libMesh/libmesh/releases/download/v1.0.0/libmesh-1.0.0.tar.gz
tar -xzvf libmesh-1.0.0.tar.gz
cd libmesh-1.0.0/ && ./configure --prefix=/usr && make && sudo make install
