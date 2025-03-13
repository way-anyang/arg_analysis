# arg_analysis

需要提供两个库：gggenes_mge_preparation.txt   &   gggenes_mge_preparation.txt
其中，两个库都必须含有contig	gene	start	end	strand 五项数据
脚本会自动检索两个库，生成gggenes的准备文件

使用：
```
Rscript arg_analysis.R NDM-1
#标黄部分为希望搜索的基因（ARG），arg_analysis.R为脚本的名称
```

运行结果如下：

<img width="338" alt="image" src="https://github.com/user-attachments/assets/6736c617-374c-468e-855f-9e8de8914966" /> <img width="294" alt="image" src="https://github.com/user-attachments/assets/d8e0fca7-bcd2-410b-b649-92778aa4b7a2" />


代码可用性：
#一般可使用，但若是含有此基因的contig数量太大，会报错，修改方法为：手动修改第190/191行的height和width
#良好的结果如左图
#由于contig数量太大，强行运行的结果如右图
