%The big idea:
%To keep the code modular, resonators and QBs are handled seperately for as long as possible.
%For example, there are seperate creation/annihilation operators for QBs and resonators.
%Later, these get combined with the Kronecker tensor product to form operators for the whole system.


%constants
global N_ph = 2; %consider up to this many photons
global N_res = 4; %number of resonators to consider
%For debugging purposes, I've put in simple numbers for the following constants
global h_bar = 1;%1.05457148e-34;
global w0 = 10;%2*pi*5e9; %resonator resonance frequency 
global wQB = 10;%2*pi*5e9; %QB resonance frequency
global wInt = 1;%2*pi*250e6; %interaction between QB and resonators
global Wcoup = 1;%2*pi*250e6; %coupling between resonators
global N_resonator_states; %run generate_resonator_states to populate
global resonator_states; %run generate_resonator_states to populate



%generates a matrix that contains the possible Fock states in rows; the number in row n is the number of photons in resonator n
%this function needs to be called before anything else
%for example, for the 2 photon, 4 resonator case, resonator_states gets populated with
%     0   0   0   0
%     1   0   0   0                                                                                                                             
%     0   1   0   0                                                                                                                             
%     0   0   1   0                                                                                                                             
%     0   0   0   1                                                                                                                             
%     2   0   0   0                                                                                                                             
%     1   1   0   0                                                                                                                             
%     1   0   1   0                                                                                                                             
%     1   0   0   1                                                                                                                             
%     0   2   0   0                                                                                                                             
%     0   1   1   0
%     0   1   0   1
%     0   0   2   0
%     0   0   1   1
%     0   0   0   2
function generate_resonator_states
  global N_ph;
  global N_res;
  global N_resonator_states;
  global resonator_states;

  N_resonator_states = 1;
  resonator_states = zeros(1,N_res);

  for ph = 1:N_ph %go over all numbers of photons except for zero photon case, which is already covered
    N_resonator_states += nchoosek(N_res+ph-1,ph);
    for st = 1:size(resonator_states)(1) %go over all the states we've already enumerated
      for res = 1:N_res %go over all the resonators
	%come up with a new state by adding a photon to the 'res'th resonator of a previously listed state
	temp = resonator_states(st,:);
	temp(res) += 1; 
	%however, this state may already have made it in to the list

	already_there = 0; %check to see if this state is already in the list
	for k = 1:size(resonator_states)(1)
	  if isequal(resonator_states(k,:),temp)
	    already_there = 1;
	  end
	end

	if already_there == 0 %add it if it's not
	  resonator_states = [resonator_states; temp];
	end
      end
    end
  end
end

% takes a population like [1,0,0,1] (one photon in first and last resonators) and returns the equivilent state: [0;0;0;0;0;0;0;0;1;0;0;0;0;0;0]
% only works for the values in the array of states, not superpositions of them
function s = density_to_resonator_state(density)
  global N_resonator_states;
  global resonator_states;
  s = zeros(N_resonator_states, 1);
  s(resonator_state_index(density)) += 1;
end

%reverse of previous function; will work on fractional states
function s = resonator_state_to_density(state)
  global N_resonator_states;
  global N_res;
  global resonator_states;
  s = zeros(1, N_res);
  
  for n = 1:N_resonator_states
    s += state(n) * resonator_states(n,:);
  end
end

function ind = resonator_state_index(state)
  global resonator_states;
  for n = 1:size(resonator_states)(1)
    if isequal(resonator_states(n,:),state)
      ind = n;
    end
  end
end


% takes a population like [1,1,0,0] (first two QBs in excited state, second two in ground state) and returns the equivilent state
% only works for the values in the array of states, not superpositions of them
function s = density_to_QB_state(density)
  global N_res; %that's how many QBs we have

  s = [1];

  for n = 1:N_res
    if density(n) == 1
      s = kron(s,[0;1]);
    else
      s = kron(s,[1;0]);
    end
  end
end

%reverse of previous function; will work on fractional states
function s = QB_state_to_density(state)
  global N_res;

end

%gives the state of the whole system from populations in the resonators and QBs
function s = density_to_state(res_density,QB_density)
  s = kron(density_to_resonator_state(res_density),density_to_QB_state(QB_density));
end

%gives the expected number of photons in each resonator for an given state
function density = state_to_resonator_desnity(state)
  global N_resonator_states;
  global N_res;
  global resonator_states;

  %the following commented out code is faster than the uncommented out code but only works in the case where all QBs are in gnd state
%    density = zeros(1,N_res);
%  
%    for n = 1:N_resonator_states;
%      prob = abs(state * density_to_state(resonator_states(n,:),[0,0,0,0]))^2; %probability we're in the 'n'th resonator fock state
%      density += resonator_states(n,:) * prob;
%    end

  %find density using expectation values of the number operators
  for n = 1:N_res
    N_op = kron(resonator_number(n),eye(2^N_res));
    density(n) = real(conj(state)*(N_op*transpose(state))); %use 'real' because small small imag part remains
  end
end

function ind = resonator_state_index(state)
  global resonator_states;
  for n = 1:size(resonator_states)(1)
    if isequal(resonator_states(n,:),state)
      ind = n;
    end
  end
end

%next two functions are helper functions; they return the new state after adding or removing a photon form a given resonator
function st = inc_resonator_state(state, res_index)
  st = state(res_index) += 1;
  st = state;
end

function st = dec_resonator_state(state, res_index)
  state(res_index) -= 1;
  st = state;
end

%HAD ANNIHILATION AND CREATION BACKWARDS; FIXED NOW, BUT FIGURE OUT WHY THIS WORKS AND CHANGE COMMENTS TO FIT
function m = resonator_annihilation(res_index)
  global N_resonator_states;
  global N_res;
  global N_ph;
  global resonator_states;

  m = zeros(N_resonator_states, N_resonator_states);

  for state = 1:(N_resonator_states-nchoosek(N_res+N_ph-1,N_ph)) %only consider the states with less than N_ph photons
    m(state,resonator_state_index(inc_resonator_state(resonator_states(state,:),res_index))) = sqrt(1+resonator_states(state,:)(res_index));
  end

end

function m = resonator_creation(res_index)
  global N_resonator_states;
  global N_res;
  global N_ph;
  global resonator_states;

  m = zeros(N_resonator_states, N_resonator_states);

  for state = 2:(N_resonator_states) %only consider the states with one or more photons
    if resonator_states(state,:)(res_index) > 0 %can't lower from a state that has nothing in it
      m(state,resonator_state_index(dec_resonator_state(resonator_states(state,:),res_index))) = sqrt(resonator_states(state,:)(res_index));
    end
  end

end

function m = resonator_number(res_index)
  m = resonator_creation(res_index)*resonator_annihilation(res_index);
end

%the creation operator for the 'QB_index'th QB
%the idea is that we just kroneker product the 2x2 creation matrix for the QB with identity QBs for the others
function m = QB_creation(QB_index)
  global N_res; %that's how many QBs we have

  m = [1];

  for n = 1:N_res
    if n == QB_index
      m = kron(m,[0,0;1,0]);
    else
      m = kron(m,eye(2));
    end
  end
end

function m = QB_annihilation(QB_index)
  global N_res; %that's how many QBs we have

  m = [1];

  for n = 1:N_res
    if n == QB_index
      m = kron(m,[0,1;0,0]);
    else
      m = kron(m,eye(2));
    end
  end
end

function m = QB_number(QB_index)
  m = QB_creation(res_index)*QB_annihilation(QB_index);
end

function m = total_QB_number
  global N_res;
  m = QB_creation(1)*QB_annihilation(1);
  for n = 2:N_res
    m += QB_creation(n)*QB_annihilation(n);
  end
end

%hamiltonian for just the resonator part of the system
function H = resonator_hamiltonian
  global N_ph;
  global N_res;
  global N_resonator_states;
  global H;
  global h_bar;
  global w0;
  global Wcoup;

  H = zeros(N_resonator_states, N_resonator_states);

  for n = 1:N_res
    H += h_bar*w0*resonator_creation(n)*resonator_annihilation(n); %yeilds the diagonal elements
    %next two 'H +=' lines add in the coupling between the resonators -- the off diagonal elements
    if n > 1
      H += h_bar*Wcoup*resonator_creation(n-1)*resonator_annihilation(n);
    end
    if n < N_res
      H += h_bar*Wcoup*resonator_creation(n+1)*resonator_annihilation(n);
    end
  end
end

%the hamiltonian for the whole system, with 'interacting' being the lists of QBs coupled to their resonators
%For example, interacting=[1,3,4] would mean QBs 1,3, and 4 are coupled to their resonators
function H = hamiltonian(interacting=[])
  global N_ph;
  global N_res;
  global h_bar;
  global w0;
  global wQB;
  global wInt;
  global Wcoup;
  global N_resonator_states;
  global resonator_states;

  H_res = resonator_hamiltonian;
  I_res = eye(N_resonator_states);
  I_qb  = eye(2^N_res);
  N_QB = total_QB_number;

  H = kron(H_res, I_qb) + kron(I_res,h_bar*wQB*N_QB);%this is the hamiltonain not including interaction between QBs and resonators
  %CHECK FOR BUG: DON'T ALLOW PHOTON TO GO IN TO QB IF QB IS ALREADY EXCITED
  for n = 1:size(interacting)(2)
    index = interacting(n);
    H += h_bar*wInt*(kron(resonator_annihilation(index), QB_creation(index)) + kron(resonator_creation(index), QB_annihilation(index)));
  end

end


function decompose(resonator_hamiltonian)
  global N

end

%evolves 'initial_state' for time 't' with indicies of interacting QBs given in 'interacting'
%NOT USEFUL AT PRESENT FOR USE AT MULTIPLE TIMES; TAKE ARGUMENTS LIKE LINSPACE???
function psi = evolve(initial_state,t,interacting=[])
  global H0;
  global Hi1;
  global Hi2;
  global h_bar;

  H = hamiltonian(interacting);

  psi = expm(-i*(H*t)/h_bar)*initial_state;
end;