require File.join(File.dirname(__FILE__), '../../tap_test_helper.rb') 
require 'ms/in_silico/fragment'

class FragmentTest < Test::Unit::TestCase
  include Ms::InSilico
  acts_as_tap_test

  #
  # headers test
  #
  
  def test_headers_returns_a_hash_of_spec_data
    frag = Fragment.new :charge => 1, :series => ['b']
    spec = Spectrum.new('TVQQEL', 'H', 'HO')
    
    headers = frag.headers(spec)
    assert_equal 1, headers[:charge]
    assert_equal 'H', headers[:nterm]
    assert_equal 'HO', headers[:cterm]
    assert_in_delta 717.377745628191, headers[:parent_ion_mass], 0.000000000001
    assert_equal ['b'], headers[:series]
  end
  
  #
  # process test
  #
  
  def test_process_returns_data_and_a_hash_of_headers
    data, headers = Fragment.new.process('TVQQEL')
    
    assert_equal %w{
      102.054954926291
      132.101905118891
      201.123368842491
      261.144498215091
      329.181946353891
      389.203075726491
      457.240523865291
      517.261653237891
      586.283116961491
      616.330067154091
      699.367180941891
      717.377745628191
    }, data.collect {|mass| mass.to_s }
    
    assert_equal 1, headers[:charge]
    assert_equal 'H', headers[:nterm]
    assert_equal 'HO', headers[:cterm]
    assert_in_delta 717.377745628191, headers[:parent_ion_mass], 0.000000000001
    assert_equal ['y', 'b'], headers[:series]
  end
end