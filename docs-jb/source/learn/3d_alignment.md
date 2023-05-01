# Structural Alignment

As mentioned earlier, to assess the structural similarity between two objects in physical space using **RMSD** 
calculation between their point clouds, the objects must first be aligned to each other.
The general problem of finding a spatial transformation (e.g. translation, rotation, scaling) to best
align two point clouds is called **point-set registration**, which has applications in various fields,
such as computer vision, pattern recognition, and robotics. In the general scenario, 
only the coordinates of the points in each point cloud is known, but not the pairwise correlation (correspondence) 
between the points of the two clouds. In such a case, the alignment algorithm must be able to find the best alignment
resulting from all possible correlations between the two sets.

However, in most CADD related scenarios, the correlation between the point clouds is known. 
The point clouds usually correspond to the coordinates of the atoms of a molecule in different conformations.
Thus, since the identity of each atom with a given coordinates is known in both conformations, the only 
valid correspondence between the two point clouds is already available. The reduced problem of aligning two
point clouds when the pairwise correlation of the two sets is known, is called **correspondence-based registration**.

```{note}
There are cases where no clear one-to-one correspondence between the points exists, for example when
comparing the sructures of two different molecules. Nevertheless, in most of such cases, methods other than
point-set registration exist that are more reliable and faster to find the suited correlation. 
After the correlation is obtained, a correspondence-based method can be used for the actual alignment. 
```

Another restriction on the problem of structural alignment of molecules is that not all spatial transformations
are allowed. Imagine a general example of point-set registration in computer vision and robotics: 
A robot is looking at an object through its cameras. In order to identify the object, it tries
to match the image of the object it captured, with labeled images in a database. To do so, it aligns the 
captured image to each of the images it has in its database, scores each alignment, and picks the label with the
highest score. In this case, it makes sense for the alignment algorithm to stretch, compress and scale the 
point cloud of the captured image, in order to account for the differences in perspectives between the captured
image of the object, and the image of that object in the database. This is because the point clouds being compared
were obtained from the projection of the object into a lower dimensional space, i.e. from the 3D physical space to 
the 2D space of the image. In contrast, the atomic coordinates data in CADD applications describe the exact 
location of each atom of the molecule in three-dimensional space. Transformations such as scaling and shearing that
do not conserve the distances between the points (here atoms) must not be allowed, since such transformations 
essentially create a different conformation. A geometric transformation that conserves the euclidean distance
between all points in space is called a **rigid transformation**, and methods to align two point clouds only 
via rigid transformations are called **rigid registration**. There is another technicality: rigid transformations
include translation, rotation and reflection (also called improper rotation). 
However, reflections should not be allowed for chiral objects, 
since it will flip the chirality of that object. Chiral objects are those that are not superposable on their mirror images. 
This includes many biological molecules, such as proteins. 

Transformations that only include translation and rotation are called **proper rigid transformations**, and the problem of
finding such a transformation that best aligns two sets of points is known as the **Wahba's problem**,{cite}`WahbaProblem` which is
essentially a restricted form of the **[orthogonal procrustes](https://docs.scipy.org/doc/scipy/reference/generated/scipy.linalg.orthogonal_procrustes.html#scipy.linalg.orthogonal_procrustes) problem** in linear algebra. 

Several solutions to the Wahba's problem exist, including those based on quaternion mathematics, such as the 
Davenport's [Q-method](https://ntrs.nasa.gov/citations/19680021122) and several other variations.{cite}`Berthold,Coutsias,Kneller,Theobald` 