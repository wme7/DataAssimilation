function [col] = nextcolor(color);
  colors = ['b' 'g' 'r' 'c' 'm' 'y' 'k'];
  
  for i = 1 : 7
    if strcmp(color, colors(i))
      break;
    end
  end
  
  col = colors(mod(i, 7) + 1);
  
  return;
