Gem::Specification.new do |s|
  s.name = "ms-in_silico"
  s.version = "0.1.0"
  s.author = "Simon Chiang"
  s.email = "simon.a.chiang@gmail.com"
  s.homepage = "http://mspire.rubyforge.org/projects/ms-in_silico/"
  s.platform = Gem::Platform::RUBY
  s.summary = "ms-in_silico task library"
  s.require_path = "lib"
  s.test_file = "test/tap_test_suite.rb"
  s.rubyforge_project = "mspire"
  s.has_rdoc = true
  s.add_dependency("tap", ">= 0.11")
  s.add_dependency("molecules", ">= 0.1.0")
  
  s.extra_rdoc_files = %W{
    README
    MIT-LICENSE
  }
  
  s.files = %W{
    lib/ms/in_silico.rb
    lib/ms/in_silico/digest.rb
    lib/ms/in_silico/digester.rb
    lib/ms/in_silico/fragment.rb
    lib/ms/in_silico/spectrum.rb
    tap.yml
  }
end