!HAMILTONIAN OF PHONONS: W0*(number of phonons)
  htmp = w0_ph*(iph - 1)
  Hv(i-MpiIshift) = Hv(i-MpiIshift) + htmp*vin(i)
