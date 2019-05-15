# Use an official Python runtime as a base image
#FROM ubuntu:17.04
FROM ubuntu:17.10

# Set the working directory to /root
WORKDIR /root

# Copy the current directory contents into the container at /app
ADD pkglst /root

# Install required packages as specified in pkglst
RUN apt-get update
RUN apt-get --assume-yes upgrade
RUN apt-get --assume-yes install $(cat pkglst)
RUN pip install jupyter
