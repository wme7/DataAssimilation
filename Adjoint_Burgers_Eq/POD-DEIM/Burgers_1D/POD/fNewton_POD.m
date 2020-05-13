function f_red=fNewton_POD(yOld_red,y_red,f_red)

global Bred Upod Cred nu dt h

y_full = Upod*y_red;

f_red = (1/dt)*y_red - (1/dt)*yOld_red + 0.5*Bred*(y_full).^2 + nu*Cred*y_red - h*f_red;

end