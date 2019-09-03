# This is a pull request, finish.
if [ "_$TRAVIS_PULL_REQUEST" != "_false" ] ;then
  echo "This is a pull request, do nothing."
  exit 0;
fi
# build doc if and only if master, develop, xxx-autodoc, and tag
feature_branch=${TRAVIS_BRANCH%-autodoc}

if [ "_$TRAVIS_BRANCH" == "_master" ]; then
  echo "This is the master branch, deploy docs."
elif [ "_$TRAVIS_BRANCH" == "_develop" ]; then
  echo "This is the develop branch, deploy docs."
elif [ "_${feature_branch}" != "_${TRAVIS_BRANCH}" ]; then
  echo "This is an auto-documented branch, deploy docs."
elif [ -n "$TRAVIS_TAG" ]; then
  echo "This is a versioned tag, deploy docs."
else
  echo "Do nothing."
  exit 0
fi


openssl aes-256-cbc -K $encrypted_0f0c7c69c924_key -iv $encrypted_0f0c7c69c924_iv -in ${ROOTDIR}/.travis_scripts/id_rsa.enc -out ~/.ssh/id_rsa -d

chmod 600 ~/.ssh/id_rsa
echo -e "Host github.com\n\tStrictHostKeyChecking no\n" >> ~/.ssh/config

git clone git@github.com:${TRAVIS_REPO_SLUG} dcore-repo
cd dcore-repo
git checkout gh-pages
if [ ${feature_branch} != ${TRAVIS_BRANCH} ]; then
  mkdir -p $feature_branch
  cp -r ${ROOTDIR}/dcore_doc/* $feature_branch
  git add $feature_branch
elif [ "_${TRAVIS_BRANCH}" == "_develop" ]; then
  mkdir -p develop
  cp -r ${ROOTDIR}/dcore_doc/* develop
  git add develop
elif [ "_${TRAVIS_BRANCH}" == "_master" ]; then
  mkdir -p master
  cp -r ${ROOTDIR}/dcore_doc/* master
  git add master
elif [ -n ${TRAVIS_TAG}]; then
  mkdir -p ${TRAVIS_TAG}
  cp -r ${ROOTDIR}/dcore_doc/* ${TRAVIS_TAG}
  git add ${TRAVIS_TAG}
else
  echo "The deploy script failed to solve where to install documents. The script has some mistake."
  echo "\$TRAVIS_BRANCH: $TRAVIS_BRANCH"
  echo "\$TRAVIS_TAG: $TRAVIS_TAG"
  echo "\$TRAVIS_PULL_REQUEST: $TRAVIS_PULL_REQUEST"
  echo "\$feature_branch: $feature_branch"
  exit 1
fi

git config --global user.email ""
git config --global user.name "DCore"
git commit -m "Update by TravisCI (\\#${TRAVIS_BUILD_NUMBER})"
ST=$?
if [ $ST == 0 ]; then
  git push origin gh-pages:gh-pages --follow-tags > /dev/null 2>&1
fi

