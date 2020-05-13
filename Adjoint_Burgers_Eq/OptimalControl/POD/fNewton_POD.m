function f_red=fNewton_POD(yOld_red,y_red,f_red, control)

global Bred Upod Cred nu dt h Mred1

y_full = Upod*y_red;

f_red = (1/dt)*y_red - (1/dt)*yOld_red + 0.5*Bred*(y_full).^2 + nu*Cred*y_red - h*f_red - Mred1*control;

end