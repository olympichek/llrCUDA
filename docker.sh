sudo docker run -v $PWD:/build -w /build -u $UID:$UID --rm --gpus=all -ti nvidia/cuda:8.0-cudnn5-devel-ubuntu16.04 "$@"
