#!/usr/bin/env bash
set -euo pipefail # Stop the script if any commands fail

##----------------------------------------------------------------------------##
# Create a basic empty project directory structure
##----------------------------------------------------------------------------##
# Version: 0.1
#
# USAGE:
#    bas-project-set-up sh <github_username> <repo_name>
#
#    NOTE: You will be asked to enter your github password
#
##----------------------------------------------------------------------------##

# get user input
MY_GIT_USERNAME=${1}
MY_PROJECT_NAME=${2}

# Make base project
echo "LOG: Start making basic project directory: ${MY_PROJECT_NAME} GutHub Username: ${MY_GIT_USERNAME}"

mkdir -v ${MY_PROJECT_NAME}

# Move into project directory
cd ${MY_PROJECT_NAME}

# Add a README.md to the project directory
echo "LOG: Adding a REAMDE.md to the new project directory"
touch README.md

# Add a basic header to the README,md
echo "# ${MY_PROJECT_NAME}" >> README.md

# Make project subdirectories
echo "LOG: Start making project subdirectories"

mkdir -v scripts analysis docs data

# Add a blank README.md to all directories
# using a for loop
for my_directory in scripts analysis docs data;do
  touch ${my_directory}/README.md
  echo "# ${my_directory}" >> ${my_directory}/README.md
done

echo "LOG: Done making basic project subdirectories: ${MY_PROJECT_NAME}"

# make public repo on github
echo "LOG: Creating private github repo"
echo "WARNING: You will be asked to enter your github password"

curl -u ${MY_GIT_USERNAME} \
https://api.github.com/user/repos \
-d "{\"name\":\"${MY_PROJECT_NAME}\",\"public\":\"false\"}"

# make it githubable
echo "LOG: Running git init"
git init

echo "LOG: Pushing the local repository to GitHub: running: git remote add origin \"https://github.com/${MY_GIT_USERNAME}/${MY_PROJECT_NAME}.git"

git remote add origin "https://github.com/${MY_GIT_USERNAME}/${MY_PROJECT_NAME}.git"

# add all and commit and puch to github
git add .
git status
git commit -m "first commit empty project"
git push --set-upstream origin master

# END
echo "LOG: Done creating new basic project directory"