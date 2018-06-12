# GA_model
R script for random network generation model incorporating the grouped attachment (GA) to resemble real-world topological properties of network motifs

## Reference
Jaejoon Choi, and Doheon Lee, "Topological motifs populate complex networks through grouped attachment", Under review at Scientific Reports

## Function Description : ga.game
Function - GA random network generation
### Input
- n - Integer. Number of nodes. (n>0)
- p - Double. Probability of edge generation. (0<p<1)
- q - Double. Groupness probability. (p<q<1)
- directed - Logical. Whether to create a directed graph. (Default : FALSE)
- revised - Logical. Whether to use revised GA model. (Default : FALSE)
- pref - Logical. Whether to use preferential GA model. (Default : FALSE)
- alpha - Initial attractiveness of the nodes. It only works if pref=TRUE. (Default : 1)
- a - Power of the preferential attachment. It only works if pref=TRUE. (Default : 1)
### Output
- An graph (igraph) object

## Contact
For any problem, please contact
- Jaejoon Choi : jjchoi@biosoft.kaist.ac.kr
- Doheon Lee : dhlee@kaist.ac.kr
