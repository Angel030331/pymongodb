# Python-SQL-Project-19-Aug-2024

## Author: Angel Wong On Ki
## Contributor: Angel Wong On Ki

This is a self-initiated project for Data Science and Machine Learning bootcamp 2024.
The aim of the project is to consolidate and use the knowledge we have learned from the Python Basics and Database Management Basics to build our own project.

My project is a database management system with MongoDB and Python interface.
The database serve serveral purposes.
 
## Usage

To use this docker image, follow these steps:

### Method 1: Clone the repository and build the docker image
After cloning the github repository, locate to the path where the Dockerfile locates.
``` bash
cd /path/to/dir
```
Run the docker compase command. Please make sure docker is installed in your system.
``` bash
docker build -t pymongodb .
docker images
```
Run the docker
``` bash
docker run -d <Image_ID>
docker ps
```
Configuration
``` bash
docker-compose -f mongodb.yml up
```

### Method 2: Pull the docker from DockerHub
