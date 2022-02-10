# Dockerfile for GitHub Actions (unit test)
FROM shinaoka/triqs3_all
COPY .github_scripts/entrypoint.sh /
COPY . /var/dcoretest
RUN ["chmod", "+x", "/entrypoint.sh"]
WORKDIR /var/dcoretest
RUN ["git", "clean", "-xdf"]
ENTRYPOINT ["/entrypoint.sh"]
#CMD ["bash"]
