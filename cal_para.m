clear all

nr = [1e2];
nr5a5 = 10;


ka = [0.1 0.2 10 0.1 0.2 10 20 0.1];
for ii = 1:length(ka)
  nr = [nr ka(ii)/nr(ii)];
  nr(ii)*nr(ii+1)
end
a = nr(2:2:end)
nr = nr(1:2:end)

ka = [nr5a5 0.1 0.2 0.4 10 0.1 0.2 10];
nr = [nr(end)];
for ii = 1:length(ka)
  nr = [nr ka(ii)/nr(ii)];
  nr(ii)*nr(ii+1)
end
a = nr(2:2:end)
nr = nr(1:2:end)
