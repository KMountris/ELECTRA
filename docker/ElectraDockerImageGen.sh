#!/bin/bash

# Read image version number from user
echo
echo \*\*\*\* ELECTRA Docker Image Generator \*\*\*\*
echo
echo Hello $USER and welcome to to the ELECTRA docker image generator.
echo Please provide the desired version number in the format major.minor.tweak \(eg., 1.0.0\)
read -p 'ELECTRA version (x.x.x): ' version
echo
echo Generating Docker image: electra-docker:$version ...

# Clean system from previous images and containers
docker system prune --force
docker volume prune --force

# Create docker image
docker image build -t electra-docker:$version .
docker save electra-docker:$version | gzip > electra-docker_$version.tar.gz

echo
echo Generated successfuly Docker image: electra-docker:$version
