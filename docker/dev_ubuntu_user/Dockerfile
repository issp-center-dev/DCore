# Use an official Python runtime as a base image
FROM shinaoka/dev_ubuntu:latest

# Become a mere user
RUN adduser testuser
USER testuser
ENV HOME=/home/testuser
WORKDIR $HOME
