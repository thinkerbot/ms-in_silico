require File.join(File.dirname(__FILE__), '../../tap_test_helper.rb') 
require 'ms/in_silico/digest'

class Ms::InSilico::DigestTest < Test::Unit::TestCase
  acts_as_script_test

  def test_digest_documentation
    script_test(File.dirname(__FILE__) +  "../../../../") do |cmd|
      cmd.check "documentation", %q{
% rap digest MIVIGRSIVHPYITNEYEPFAAEKQQILSIMAG --+ dump --no-audit
  I[:...:]             digest MIVIGRSIVHP... to 3 peptides
# date: :...:
--- 
ms/in_silico/digest (:...:): 
- - MIVIGR
  - SIVHPYITNEYEPFAAEK
  - QQILSIMAG
}
    end
  end
end