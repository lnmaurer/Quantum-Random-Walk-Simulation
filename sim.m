%constants
global N_ph = 2; %consider up to this many photons
global N_res = 4; %number of resonators to consider
global h_bar = 1.05457148e-34;
global w0 = 2*pi*5e9; 
global Wcoup = 2*pi*250e6;
%  global N_states = 15;
%  global states = [0,0,0,0;
%  		 1,0,0,0;
%  		 0,1,0,0;
%  		 0,0,1,0;
%  		 0,0,0,1;
%  		 2,0,0,0;
%  		 1,1,0,0;
%  		 1,0,1,0;
%  		 1,0,0,1;
%  		 0,2,0,0;
%  		 0,1,1,0;
%  		 0,1,0,1;
%  		 0,0,2,0;
%  		 0,0,1,1;
%  		 0,0,0,2;];


%generates a matrix that contains the possible Fock states in rows; the number in row n is the number of photons in resonator n
function states = generate_states
  global N_ph;
  global N_res;
  global N_states;
  global states;

  N_states = 1;
  states = [0,0,0,0];

  for ph = 1:N_ph %go over all numbers of photons except for zero photon case, which is already covered
    N_states += nchoosek(N_res+ph-1,ph);
    for st = 1:size(states)(1) %go over all the states we've already enumerated
      for res = 1:N_res %go over all the resonators
	%come up with a new state by adding a photon to the 'res'th resonator of a previously listed state
	temp = states(st,:);
	temp(res) += 1; 
	%however, this state may already have made it in to the list

	already_there = 0; %check to see if this state is already in the list
	for k = 1:size(states)(1)
	  if isequal(states(k,:),temp)
	    already_there = 1;
	  end
	end

	if already_there == 0 %add it if it's not
	  states = [states; temp];
	end
      end
    end
  end
end

% takes a density like [1,0,0,1] and returns the equivilent state: [0;0;0;0;0;0;0;0;1;0;0;0;0;0;0]
% only works for the values in the array of states, not for fractional states
function s = density_to_state(density)
  global N_states;
  global states;
  s = zeros(N_states, 1);
  s(state_index(density)) += 1;
end

%reverse of previous function; will work on fractional states
function s = state_to_density(state)
  global N_states;
  global N_res;
  global states;
  s = zeros(1, N_res);
  
  for n = 1:N_states
    s += state(n) * states(n,:);
  end
end

function ind = state_index(state)
  global states;
  for n = 1:size(states)(1)
    if isequal(states(n,:),state)
      ind = n;
    end
  end
end

function st = inc_state(state, res_index)
  st = state(res_index) += 1;
  st = state;
end

function st = dec_state(state, res_index)
  state(res_index) -= 1;
  st = state;
end

function m = creation(res_index)
  global N_states;
  global N_res;
  global N_ph;
  global states;

  m = zeros(N_states, N_states);

  for state = 1:(N_states-nchoosek(N_res+N_ph-1,N_ph)) %only consider the states with less than N_ph photons
    m(state,state_index(inc_state(states(state,:),res_index))) = sqrt(1+states(state,:)(res_index));
  end

end

function m = annihilation(res_index)
  global N_states;
  global N_res;
  global N_ph;
  global states;

  m = zeros(N_states, N_states);

  for state = 2:(N_states) %only consider the states with one or more photons
    if states(state,:)(res_index) > 0 %can't lower from a state that has nothing in it
      m(state,state_index(dec_state(states(state,:),res_index))) = sqrt(states(state,:)(res_index));
    end
  end

end

function H = hamiltonian
  global N_ph;
  global N_res;
  global N_states;
  global H;
  global h_bar;
  global w0;
  global Wcoup;

  H = zeros(N_states, N_states);

  for n = 1:N_res
    H += h_bar*w0*annihilation(n)*creation(n);
    if n > 1
      H += h_bar*Wcoup*annihilation(n-1)*creation(n);
    end
    if n < N_res
      H += h_bar*Wcoup*annihilation(n+1)*creation(n);
    end
  end
end

function decompose(hamiltonian)
  global N

end

function rho = density(initial_state,t)
  global H0;
  global Hi1;
  global Hi2;
  global h_bar;

  H = hamiltonian;

  rho = expm(-i*(H*t)/h_bar)*initial_state*expm(i*(H*t)/h_bar);
end;