#!/bin/bash

if [ ! `which gh` ]
then
    echo "[ERROR] INSTALL github-cli: https://github.com/cli/cli/blob/trunk/docs/install_linux.md"
    exit 1
fi

if ! (($#)); then
  echo "[ERROR] {NEW_REPO_NAME} should be given!" >&2
  cat << EOF

CREATEPROJECT.SH creates a project from a template (PROJTEMP) in Github. It requires Github command line interface "GH".
USAGE:    ./createprojet.sh {NEW_REPO_NAME}
EXAMPLE:  ./createprojet.sh my_new_awesome_project
(cc) 2024, dr.seunggoo.kim@gmail.com
EOF
  exit 2
fi


ProjName=${1}
gh repo create ${ProjName} --template projtemp --private
gh repo clone ${ProjName}
cd ${ProjName}
cat << EOF > README.md
# ${ProjName}

## directory structure

### Synchronized
\`./code\`  for functions, scripts, small (<50 MB) packages for each language

\`./docs\`  for documents (notes, presentation slides, manuscripts, other human-readable files)

\`./figs\`  to dump figure files in

\`./meta\`  for small (<50 MB) meta-data files

### Not synchronized
\`./data\`  for large (>50 MB) files. \
The folder name "data" is set NOT TO SYNCHRONIZED (see \`.gitignore\`)
EOF
rm create_project.sh
rm pull_everything.sh
mkdir data
git add .
git commit -m "Initial commit"
git push
