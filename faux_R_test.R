
xvar = c(1,2,3,4,2,1,2,4,3)
yvar = c(6,5,7,8,9,7,6,5,3)
write(t.test(xvar, yvar)$p.value, file="fake.t.test.txt")

