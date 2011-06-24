%The big idea:
%To keep the code modular, resonators and QBs are handled seperately for as long as possible.
%For example, there are seperate creation/annihilation operators for QBs and resonators.
%Later, these get combined with the Kronecker tensor product to form operators for the whole system.

%NOTES:
%-states are column vectors
%-densities/populations/whatever-you-want-to-call-them and straight up fock states are row vectors

%constants
global N_ph = 2; %consider up to this many photons
global N_res = 2; %number of resonators to consider
%For debugging purposes, I've put in simple numbers for the following constants
global h_bar = 1.05457148e-34;
global w0 = 2*pi*5e9; %resonator resonance frequency 
global wQB = 2*pi*5e9; %QB resonance frequency
global wInt = 2*pi*1/10e-9/4; %interaction between QB and resonators
global Wcoup = 2*pi*1/100e-9/4; %coupling between resonators
global N_resonator_states; %run generate_resonator_states to populate
global resonator_states; %run generate_resonator_states to populate


function ev = expectation_value(state, operator)
  ev = real(state'*(operator*state)); %use 'real' because small small imag part remains
end


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

% takes a fock state like [1,0,0,1] (one photon in first and last resonators) and returns the equivilent state: [0;0;0;0;0;0;0;0;1;0;0;0;0;0;0]
% only works for the values in the array of states, not superpositions of them
function s = fock_state_to_resonator_state(fock_state)
  global N_resonator_states;
  global resonator_states;
  s = zeros(N_resonator_states, 1);
  s(resonator_state_index(fock_state)) += 1;
end

%reverse of previous function; will work on superpositions of states
function density = resonator_state_to_density(state)
  global N_res;
  for n = 1:N_res
    density(n) = expectation_value(state, resonator_number(n));
  end
end

%returns the row number for 'fock_state' in 'resonator_states' matrix
function ind = resonator_state_index(fock_state)
  global resonator_states;
  for n = 1:size(resonator_states)(1)
    if isequal(resonator_states(n,:),fock_state)
      ind = n;
    end
  end
end


% takes a fock_state like [1,1,0,0] (first two QBs in excited state, second two in ground state) and returns the equivilent state
% doesn't work on superpositions of fock_states not superpositions of them
function s = fock_state_to_QB_state(fock_state)
  global N_res; %that's how many QBs we have

  s = [1];

  for n = 1:N_res
    if fock_state(n) == 1
      s = kron(s,[0;1]); %'n'th resonator in excited state
    else
      s = kron(s,[1;0]); %'n'th resonator in gnd state
    end
  end
end

%reverse of previous function; will work on fractional states
function s = QB_state_to_density(state)
  global N_res;
  for n = 1:N_res
    density(n) = expectation_value(state, QB_number(n));
  end
end

%gives the state of the whole system from fock_states for resonators and QBs
function s = density_to_state(res_fock,QB_fock)
  s = kron(fock_state_to_resonator_state(res_fock),fock_state_to_QB_state(QB_fock));
end

%gives the expected number of photons in each resonator for an state or states
function densities = state_to_resonator_desnity(states)
  global N_res;

  for n = 1:N_res %generate all the number operators we'll need
    N_op{n} = kron(resonator_number(n),eye(2^N_res));
  end

  densities = [];
  for state = states %iterate over the states
    for n = 1:N_res %find the denisty for the state we're on
      density(n) = expectation_value(state, N_op{n});
    end
    densities = [densities; density];%add the denisty we found to the matrix of all densities
  end
end

%gives the expected number of photons in each resonator for an given state
function densities = state_to_QB_desnity(states)
  global N_resonator_states;
  global N_res;

  for n = 1:N_res
    N_op{n} = kron(eye(N_resonator_states),QB_number(n));
  end

  densities = [];
  for state = states %iterate over the states
    for n = 1:N_res %find the denisty for the state we're on
      density(n) = expectation_value(state, N_op{n});
    end
    densities = [densities; density];%add the denisty we found to the matrix of all densities
  end
end

%next two functions are helper functions; they return the fock state resulting from adding or removing a photon form a given resonator
function st = inc_resonator_fock_state(state, res_index)
  state(res_index) += 1;
  st = state;
end

function st = dec_resonator_fock_state(state, res_index)
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
    m(state,resonator_state_index(inc_resonator_fock_state(resonator_states(state,:),res_index))) = sqrt(1+resonator_states(state,:)(res_index));
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
      m(state,resonator_state_index(dec_resonator_fock_state(resonator_states(state,:),res_index))) = sqrt(resonator_states(state,:)(res_index));
    end
  end

end

%number operator on 'res_index'th resonator
function m = resonator_number(res_index)
  m = resonator_creation(res_index)*resonator_annihilation(res_index);
end

%the creation operator for the 'QB_index'th QB
%the idea is that we just kroneker product the 2x2 creation matrix for the 'QB_index'th QB with identity QBs for the others
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

%number operator on 'QB_index'th QB
function m = QB_number(QB_index)
  m = QB_creation(QB_index)*QB_annihilation(QB_index);
end

%sum of individual QB_number operators
function m = total_QB_number
  global N_res;
  m = 0;

  for n = 1:N_res
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

%hamiltonian for just the QB part of the system
function H = QB_hamiltonian
  global h_bar;
  global wQB;

  H = h_bar*wQB*total_QB_number;
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
  H_QB  = QB_hamiltonian;
  I_res = eye(N_resonator_states);
  I_qb  = eye(2^N_res);

  H = kron(H_res, I_qb) + kron(I_res,H_QB);%this is the hamiltonain not including interaction between QBs and resonators

  %now, add in the requested interaction(s)
  for n = 1:size(interacting)(2)
    index = interacting(n);
    H += h_bar*wInt*(kron(resonator_annihilation(index), QB_creation(index)) + kron(resonator_creation(index), QB_annihilation(index)));
  end
end

%evolves 'initial_state' for time 't' with indicies of interacting QBs given in 'interacting'
function psi = evolve(state,t,interacting=[])
  global H0;
  global Hi1;
  global Hi2;
  global h_bar;

  H = hamiltonian(interacting);

  %method 1, decompose state in to the eigenvectors of 'H', then evolve those -- this seems to be fastest, especially for more than one time value
  [eig_vectors,D] = eig(H);
  eig_values = diag(D);
  psi = zeros(size(H)(1),size(t)(2)); %number of rows equal to number of states, number of columns equal to number of times
  for n = 1:size(eig_vectors)(1)
    projection = eig_vectors(:,n) * (eig_vectors(:,n)'*state); %projection of 'state' on to 'n'th eigenvector
    for m = 1:size(t)(2)
      psi(:,m) += exp(-i*eig_values(n)*t(m)/h_bar)*projection; %evolve weighted 'n'th eigenvector and add to results
    end
  end

  %method 2, exponentiate the matrix -- this seems to be slower; would need to rewrite to handle 't' being an array
%    psi = expm(-i*(H*t)/h_bar)*state;
end;

%evolves 'initial_state' to times [0,stop_time/(npoints-1),2*stop_time/(npoints-1),...,stop_time] with indicies of
%interacting QBs given in 'interacting' this function works by exponentiating a matrix; that's a high fixed time cost,
%but it only has to be done once because the times are easilly space, so this can be faster than the previous
%function if you're dealing with many evenly spaced time steps
function psi = evolve_linspaced_t(state,stop_time,npoints,interacting=[])
  global H0;
  global Hi1;
  global Hi2;
  global h_bar;

  H = hamiltonian(interacting);
  dt = stop_time/(npoints-1);

  ev = expm(-i*H*dt/h_bar);

  psi(:,1) = state;

  for n = 2:npoints
    psi(:,n) = ev*psi(:,n-1);
  end
end;

function js = HM(wait_time)
  global h_bar;
  global wInt;

  %TODO: MAKE ALL THIS GENERAL
  psi = density_to_state([0,0],[1,1]);
  psi = expm(-i*hamiltonian([1,2])*10e-9/h_bar)*psi; %swap from QBs to resonators
  psi = expm(-i*hamiltonian()*wait_time/h_bar)*psi; %wait with photons in resonators
  psi = expm(-i*hamiltonian([1,2])*10e-9/h_bar)*psi; %swap back from resonators to QBs

  p11 = abs(psi'*density_to_state([0,0],[1,1]))^2;
  p01 = abs(psi'*density_to_state([0,0],[0,1]))^2;
  p10 = abs(psi'*density_to_state([0,0],[1,0]))^2;

  js = p11/p01/p10;
end
  
function p = ps(wait_time)
  global h_bar;
  global wInt;

  %TODO: MAKE ALL THIS GENERAL
  psi = density_to_state([0,0],[1,1]);
  psi = expm(-i*hamiltonian([1,2])*10e-9/h_bar)*psi; %swap from QBs to resonators
  psi = expm(-i*hamiltonian()*wait_time/h_bar)*psi; %wait with photons in resonators
  psi = expm(-i*hamiltonian([1,2])*10e-9/sqrt(2)/h_bar)*psi; %swap back from resonators to QBs

  p(1) = abs(psi'*density_to_state([0,0],[1,1]))^2;
  p(2) = abs(psi'*density_to_state([1,0],[0,1]))^2 + abs(psi'*density_to_state([0,1],[0,1]))^2;
  p(3) = abs(psi'*density_to_state([1,0],[1,0]))^2 + abs(psi'*density_to_state([0,1],[1,0]))^2;
  p(4) = abs(psi'*density_to_state([1,1],[0,0]))^2;
  p(5) = abs(psi'*density_to_state([2,0],[0,0]))^2;
  p(6) = abs(psi'*density_to_state([0,2],[0,0]))^2;
end