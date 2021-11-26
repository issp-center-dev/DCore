FROM shinaoka/triqs3_minimum
COPY ./docker/entrypoint.sh /
RUN ["chmod", "+x", "/entrypoint.sh"]
ENTRYPOINT ["/entrypoint.sh"]
CMD ["bash"]
