# You must first have docker installed on your system.

# Clone the gemchem repository from Github
git clone https://github.com/metno/genchem.git

# Go to the directory and build the docker image from the Dockerfile
docker build -t genchem_image .

# See the image listed
docker image ls

# Use the image to run a container
# -it. The session is interactive so the user can interact and use the container
# -v.  Mounts a file space on the machine to pass files to and from the container
docker run --name genchem_container -it -v C:\directory\on\disk\:/Code/genchem/ genchem_image

# In a separate terminal you can see your docker image
docker container ls

# to remove the container
docker rm docker_container

# to remove the image
docker rmi docker_image
