## rewrite TODO sections for your project

# This is a pull request, finish.
if [ "_$TRAVIS_PULL_REQUEST" != "_false" ] ;then exit 0; fi
# This is neither master nor tag, finish.
feature_branch=${TRAVIS_BRANCH%-autodoc}
if [ "_$TRAVIS_BRANCH" != "_master" ] && [ ${feature_branch} == ${TRAVIS_BRANCH} ] && [ -z "$TRAVIS_TAG" ] ; then exit 0; fi


openssl aes-256-cbc -K $encrypted_0f0c7c69c924_key -iv $encrypted_0f0c7c69c924_iv -in ${ROOTDIR}/.travis_scripts/id_rsa.enc -out ~/.ssh/id_rsa -d
openssl aes-256-cbc -K "$encrypted_aa0e0f6aad31_key" -iv "$encrypted_aa0e0f6aad31_iv" -in ${ROOTDIR}/.travis_scripts/id_rsa.enc -out ~/.ssh/id_rsa -d

chmod 600 ~/.ssh/id_rsa
echo -e "Host github.com\n\tStrictHostKeyChecking no\n" >> ~/.ssh/config

git clone git@github.com:${TRAVIS_REPO_SLUG} dcore-repo
cd dcore-repo
git checkout gh-pages
if [ ${feature_branch} != ${TRAVIS_BRANCH} ]; then
  cp -r ${ROOTDIR}/dcore_doc $feature_branch
  git add $feature_branch
elif [ -n ${TRAVIS_TAG}]; then
  cp -r ${ROOTDIR}/dcore_doc ${TRAVIS_TAG}
  git add ${TRAVIS_TAG}
else
  cp -r ${ROOTDIR}/dcore_doc/* .
  git add .
fi

git config --global user.email ""
git config --global user.name "DCore"
git commit -m "Update by TravisCI (\\#${TRAVIS_BUILD_NUMBER})"
ST=$?
if [ $ST == 0 ]; then
  git push origin gh-pages:gh-pages --follow-tags > /dev/null 2>&1
fi

