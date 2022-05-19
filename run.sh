g++ \
-I /opt/X11/include -L /opt/X11/lib -lX11 \
*.cpp -o main \
&& rm -f \
data/phi/*.vtk \
figures/phi/*.png \
&& ./main\
&& rm -f main