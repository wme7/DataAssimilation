function obloc=obs_location(nx,nob,obs_location_type)
for p=1:nob
  if(strcmp(obs_location_type,'random'))
    obloc(p)=rand*(nx-1)+1;
  end
  if(strcmp(obs_location_type,'regular'))
    obloc(p)=p*nx/nob;
  end
  if(strcmp(obs_location_type,'regular1'))
    obloc(p)=(p-1)*nx/nob+1;
  end
  if(strcmp(obs_location_type,'partial'))
    obloc(p)=(nx-nob)/2+p-1;
  end
end
end
