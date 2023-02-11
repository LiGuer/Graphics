* Ear Clipping 
  - Purpose
    The ear clipping algorithm is a simple and efficient method for triangulating simple polygon. THe algorithm takes as input an $n$-sided polygon and produces a list of $n-2$ triangular faces that form the triangulated polygon.

  - Idea
    Identify and remove "ears" from the polygon. An ears is a triangle form by three consecutive vertices of the polygon, such that the triangle is entirely contained whitin the polygon and no other vertices of the polygon are inside the triangle.

  - Process
     
    - judge a vertice whether is a reflex or convex vertice.
      $$det(b - a, c - b) = (b_x - a_x) (c_y - b_y) - (c_x - b_x) (b_y - a_y) > 0$$

      Due to an input polygon is in most cases a list of consecutive vertices, and represent the polygon in counter clockwise order. Therefore, when walking along the boundary of the polygon the interior of it should be to the left of evert edge traversed.

      The determinate gives us as result greater than zero if the vertice $a, b, c$ form a convex angle at $b$; and less than zero otherwise.
      
      (reference: ```https://stackoverflow.com/questions/40410743/polygon-triangulation-reflex-vertex```)