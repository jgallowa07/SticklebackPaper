'''
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import numpy as np



fig, axs = plt.subplots(3, 1)

img1 = mpimg.imread("/Users/jaredgalloway/Desktop/HAPSREAL/5e_5p1.png")
img2 = mpimg.imread("/Users/jaredgalloway/Desktop/HAPSREAL/5e_5p2.png")
img3 = mpimg.imread("/Users/jaredgalloway/Desktop/HAPSREAL/5e_5p3.png")

axs[0] = plt.imshow(img1)
axs[1] = plt.imshow(img2)
axs[2] = plt.imshow(img3)

plt.show()
'''



import numpy as np
import matplotlib.pyplot as plt

w=10
h=10
fig=plt.figure(figsize=(8, 8))
columns = 4
rows = 5
for i in range(1, columns*rows +1):
    img = np.random.randint(10, size=(h,w))
    fig.add_subplot(rows, columns, i)
    plt.imshow(img)
plt.show()


