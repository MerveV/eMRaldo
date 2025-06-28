### **Node Size:**

You can choose to adjust the size of nodes based on either FDR values or the number of overlapped genes (6th column in the table). When the **FDR value is smaller**, the corresponding nodes will appear **larger**, indicating higher significance in the dataset.


### **Node Color:**

Nodes can be colored based on z-score of each selected dataset or differences of datasets. z-score is calculated as the mean of the logFC values of the genes ('hits' column in table) in selected dataset divided by the standard deviation of those values. If 'Difference' is selected, it first takes the differences of logFC values between two selected datasets, then calculates the z-score.

If the calculated **z-score is greater than 0**, the nodes in the network are colored as **red**, indicating upregulation, otherwise as **blue** representing downregulation.

### **Node Border:**

The border of nodes is drawn based on the absolute value of the **z-score multiplied by 10**.

Upon clicking or hovering on any node, it shows z-score and genes hit by that TF. 

### **Node Shape:**

Wilcoxon test is used to find if the average of logFC values of the genes hit by that TF is significantly different from 0. If it is **significanly different**, the node is drawn as **diamond**, otherwise as circular. 

### **Edge Width:**

The width of edges in the network is determined by the number of shared genes between two transcription factors (TFs). **Thicker edges** indicate a **higher number of shared genes**.

### **Edge Color:**

Mean of logFC values of genes shared by two TFs are used to color edges. If it is **greater than 0**, the edge is colored as **red**, indicating upregulation, otherwise as **blue** representing downregulation.

Upon clicking or hovering on any edge, it shows mean logFC value and shared genes between two TFs. 

### **Edge Type:**

Wilcoxon test is used to find if the average of logFC values of the shared genes by 2 TFs is significantly different from 0. If it is **significanly different**, the edge is drawn as **solid**, otherwise as **dashed**. 

The **mode of regulation** is represented with an arrow if it is an activator, and with an inhibitor if it is a repressor. Activators and inhibitors are also labeled **A** and **R**, respectively.



## **To cut the edges:**

* Jaccard distance: Represents the similarity between the common genes of two transcription factors (TFs), where 1 indicates the highest similarity and 0 indicates dissimilarity. 
* Number of shared genes: Shows the count of shared/common genes between the two TFs.
* Percentage of shared genes: Calculated as **100 x (the intersection of TFs / their union)**
* Edge's logFC: Represents the cutoff for the absolute value of the mean logFC values of shared genes between two TFs. 
* Significance: Displays the result of the Wilcoxon test. If **'Significant'** is selected, only significant edges will be shown, while **'Not Significant'** will display non-significant edges as dashed lines.

<em>Edges having greater value than selected cutoff will be drawn. </em>
