<script type="text/javascript" src="http://cdn.mathjax.org/mathjax/latest/MathJax.js?config=default"></script>

# Matrix Map
$$u_{i,j} : j*xSize+i$$

## Formula : Down side border(j = 0)
$$u_{i,j} = 0$$

## Formula : Up side border(j = 1)
$$u_{i,j} = 0$$

## Formula : Left side border(i = 0)
$$u_{i,j} = \frac{sin({\pi}*y)}{9+{\pi}^2}$$

## Formula : Up side border(i = \\(\pi\\))
$$u_{i,j} = -\frac{sin({\pi}*y)}{9+{\pi}^2}$$

### Equation : 五点差分格式
$$u_{ij}-\frac{1}{4}u_{i-1,j}-\frac{1}{4}u_{i+1,j}-\frac{1}{4}u_{i,j-1}-\frac{1}{4}u_{i,j+1}=\frac{h^2}{4}f_{i,j}$$
#### Replace
$$f_{i,j}=-cos(3x_i)sin({\pi}y_j)$$
$$x_i=h\*i$$
$$y_j=h\*j$$

<meta http-equiv="refresh" content="30">