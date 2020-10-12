# CG-Project2
圖導 project2

2020/10/10

找出projection matrix 以及 modelview matrix
但是分別有一點小問題，不影響顯示

1. projection matrix傳入之後再取出來，顯示的值不相同
projection_matrix[11] 與 projection_matrix[14] 是相反的值，剩下都一樣

2. modelview得出的矩陣，除了w,也就是modelview[15]以外，其他值都與glulookat得出的值正負號相反
我直接把他們乘上-1來解決
