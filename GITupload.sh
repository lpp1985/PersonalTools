chmod -R 755 ./*.*
#dos2unix *.*
#GIT_TRACE=1
git add . -A
git commit -m 123
git push -f origin master:master
