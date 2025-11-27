#!/bin/bash
# Script to push TE_overdispersion to GitHub
# Make sure you've created the repository on GitHub first at: https://github.com/new

git remote add origin https://github.com/Ade-omole/TE_overdispersion.git
git branch -M main
git push -u origin main

echo "Repository pushed to GitHub successfully!"

