function f_red=fNewton_PODDEIM(yOld_red,y_red,f_red, control)

global Bred Fred Cred nu dt h Mred1

f_red = (1/dt)*y_red - (1/dt)*yOld_red + 0.5*Bred*(Fred*y_red).^2 + nu*Cred*y_red - h*f_red - Mred1*control;

end