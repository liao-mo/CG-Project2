# CG-Project2
圖導 project2

2020/10/10
# 找出projection matrix 以及 modelview matrix
但是分別有一點小問題，不影響顯示

1. projection matrix傳入之後再取出來，顯示的值不相同
projection_matrix[11] 與 projection_matrix[14] 是相反的值，剩下都一樣

2. modelview得出的矩陣，除了w,也就是modelview[15]以外，其他值都與glulookat得出的值正負號相反
我直接把他們乘上-1來解決

2020/10/21
# 完成clipping，可以順利的用glVertex2f來畫出所有的牆

clipping的步驟
1. 先把所有edge從world space轉到view space，得到所有在camera前的所有edges座標
2. 建構裁切平面(四條線)的方程式，這邊是以x為右邊，z為前面的座標系統來實作
3. 因為view space中edges的座標，是以-z作為camera前面的格式，所以要先把-z轉成z，做完clip，再把負號放回來
4. 用四個clipping plane對所有edges做clipping，沒有在可視範圍內的點就全部設為(-1,-1,-1,-1)
5. 把所有處理完的edges存入一個vector中，再交給draw wall去處理
