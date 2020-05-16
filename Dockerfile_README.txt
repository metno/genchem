# If you run on a system that is not compatible with genchem e.g. on a windows
# machine, you might consider using docker as a way of accessing genchem.
# Docker is a virtualisation software that allows the running of programmes
# across different platforms. If you have docker installed you can simply
# dowload the Dockerfile found in this repository to build a docker image and
# then from that image make a container (similar to virtual environment) in which
# to run genchem.

# You must first have docker installed on your system. You can do that by installing
# docker desktop here https://www.docker.com/products/docker-desktop. which you can
# use from the command prompt / terminal.

# Clone the gemchem repository from Github
git clone https://github.com/metno/genchem.git

# Go to the directory and build the docker image from the Dockerfile
docker build -t genchem_image .

# See the image listed
docker image ls

# Use the image to run a container
# -it. The session is interactive so the user can interact and use the container
# --name <arg>. Give the container a name
docker run --name genchem_container -it genchem_image

# In a separate terminal you can see your docker image
docker container ls

# to remove the container
docker rm docker_container

# to remove the image
docker rmi docker_image
