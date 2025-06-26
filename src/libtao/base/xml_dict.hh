#ifndef tao_base_xml_dict_hh
#define tao_base_xml_dict_hh

#include <list>
#include <iostream>
#include <typeinfo>
#include <sstream>
#include <string>
#include <vector>
#include <boost/optional.hpp>
#include <pugixml.hpp>
#include <libhpc/debug/assertions.hh>

namespace tao {
   using namespace pugi;

   class bad_option
      : public hpc::exception
   {
   };

   ///
   /// Exception thrown when there is any trouble parsing XML.
   ///
   class bad_xml
      : public bad_option
   {
   };

   ///
   /// A direct XML based options dictionary. The standard dictionary
   /// uses option predeclaration to allow program help and unified
   /// defaults. However, sometimes it is more convenient to not have
   /// to require that predeclaration. xml_dict provides this.
   ///
   class xml_dict
   {
   public:

      ///
      /// Default constructor.
      ///
      xml_dict();

      ///
      /// XML base constructor.
      ///
      xml_dict( xml_node root );

      xml_dict( std::istream& strm,
		const std::string& xpath_root = std::string() );

      xml_dict( const xml_dict& src );

      ///
      /// Destructor.
      ///
      ~xml_dict();

      ///
      /// Read an XML from file. If called multiple times, each subsequent
      /// XML will be merged with the existing dictionary.
      ///
      /// @param[in] filename The filename of the XML file to load.
      /// @param[in] xpath_root Optional XPath root to use.
      /// @throws bad_xml If any error occurs during XML parsing.
      ///
      void
      read( const std::string& filename,
	    const std::string& xpath_root = std::string() );

      ///
      /// Read an XML from a C++ stream. If called multiple times, each
      /// subsequent XML will be merged with the existing dictionary.
      ///
      /// @param[in] stream The C++ stream to read from.
      /// @param[in] xpath_root Optional XPath root to use.
      /// @param[in] filename Optional filename, used for error reporting.
      /// @throws bad_xml If any error occurs during XML parsing.
      ///
      void
      read( std::istream& stream,
	    const std::string& xpath_root = std::string(),
	    const std::string& filename = std::string() );

      ///
      /// Test for existence of an option.
      ///
      /// @param[in] path The option path to test for.
      /// @returns Boolean describing the existence of the option.
      ///
      bool
      has( const std::string& path ) const;

      ///
      /// Extract and return an option value.
      ///
      /// @tparam T The type of the option.
      /// @param[in] path The option path to get the value of.
      /// @returns The value of the option.
      /// @throws bad_option Thrown when the option path does not exist.
      ///
      template< class T >
      T
      get( const std::string& path ) const
      {
	 return _coerce<T>( _get_node( path ).first_child().value() );
      }

      ///
      /// Extract and return an option value in an optional
      /// structure.
      ///
      /// @tparam T The type of the option.
      /// @param[in] path The option path to get the value of.
      /// @returns The optional value of the option, depending on whether
      ///          the option name exists.
      ///
      template< class T >
      boost::optional<T>
      opt( const std::string& path ) const
      {
	 xml_node node = _get_node( path, false );
	 if( node )
	    return _coerce<T>( node.first_child().value() );
	 else
	    return boost::none;
      }

      ///
      /// Extract and return an option value. If the option does not
      /// exist default_value is returned.
      ///
      /// @tparam T The type of the option.
      /// @param[in] path The option path to get the value of.
      /// @param[in] default_value The value to return if the option path
      ///                          cannot be found.
      /// @returns The value of the option.
      ///
      template< class T >
      T
      get( const std::string& path,
	   const T& default_value ) const
      {
	 xml_node node = _get_node( path, false );
	 if( node )
	    return _coerce<T>( node.first_child().value() );
	 else
	    return default_value;
      }

      ///
      /// Extract and return an option value as a list. An option list
      /// in XML is interpreted as being a set of subelements. Each sub-
      /// element with a pcdata child has that pcdata value appended to
      /// the list.
      ///
      /// @tparam T The type of the elements of the list.
      /// @param[in] path The option path to get the value of.
      /// @returns The list of values.
      ///
      template< class T >
      std::list<T>
      get_list( const std::string& path ) const
      {
	 xml_node node = _get_node( path );
	 std::list<T> val;
	 for( xml_node_iterator it = node.begin(); it != node.end(); ++it )
	 {
	    if( it->first_child() && it->first_child().type() == node_pcdata )
	       val.push_back( _coerce<T>( it->first_child().value() ) );
	 }
	 return val;
      }

      ///
      /// Extract and return an option value as a list. Same as
      /// the former get_list but allows for a default value if
      /// the list option name cannot be found.
      ///
      /// @tparam T The type of the elements of the list.
      /// @param[in] path The option path to get the value of.
      /// @param[in] default_value The default value.
      /// @returns The list of values.
      ///
      template< class T >
      std::list<T>
      get_list( const std::string& path,
		const std::list<T>& default_value )
      {
	 xml_node node = _get_node( path, false );
	 if( node )
	 {
	    std::list<T> val;
	    for( xml_node_iterator it = node.begin(); it != node.end(); ++it )
	    {
	       if( it->first_child() && it->first_child().type() == node_pcdata )
		  val.push_back( _coerce<T>( it->first_child().value() ) );
	    }
	    return val;
	 }
	 else
	    return default_value;
      }

      ///
      /// Extract and return an option value as a list of attributes.
      ///
      /// @tparam T The type of the elements of the list attribute.
      /// @param[in] path The option path to get the value of.
      /// @param[in] attribute The name of the attribute to return.
      /// @returns A list of optional values. For each entry of the list
      ///          that contains an attribute accordingly named a non
      ///          'none' entry in the list will be given.
      ///
      template< class T >
      std::list<boost::optional<T> >
      get_list_attributes( const std::string& path,
			   const std::string& attribute ) const
      {
	 xml_node node = _get_node( path );
	 std::list<boost::optional<T> > val;
	 for( xml_node_iterator it = node.begin(); it != node.end(); ++it )
	 {
	    if( it->first_child() && it->first_child().type() == node_pcdata )
	    {
	       xml_attribute attr = it->attribute( attribute.c_str() );
	       if( !attr.empty() )
		  val.push_back( _coerce<T>( attr.value() ) );
	       else
		  val.push_back( boost::none );
	    }
	 }
	 return val;
      }

      xpath_node_set
      get_nodes( const std::string& xpath ) const;

      xml_node
	  get_root() const;

   protected:

      xml_node
      _find_root( xml_node& node,
		  const std::string& xpath_root ) const;

      void
      _merge( std::istream& stream,
	      const std::string& path,
	      const std::string& filename = std::string() );

      void
      _merge_node( xml_node merge_into,
		   xml_node merge_from );

      xml_node
      _get_node( const std::string& path,
		 bool except = true ) const;

      std::string
      _xform_path( const std::string& path ) const;

      template< class T >
      T
      _coerce( const std::string& value ) const
      {
	 std::stringstream ss( value );
	 T val;
	 ss >> val;
	 if( ss.fail() )
	    throw std::bad_cast();
	 return val;
      }

   protected:

      std::string _sep;
      xml_document _doc;
      xml_node _root;
   };

   template<>
   std::string
   xml_dict::_coerce( const std::string& value ) const;

}

class cli_dict
   {
   public:
      cli_dict()
      {
      }
      // User options through set trough cli options
      std::string _dataset;
      double _decmin = -90.0;
      double _decmax = 90.0;
      double _ramin  = 0.0;
      double _ramax  = 360.0;
      double _zmin   = 0.0;
      double _zmax   = 10.0;
      std::string _outfile = "output.hdf5";
      std::vector<std::string> _output_fields;

      std::string _filter_field = "stellarmass";
      std::string _filter_min = "0";
      std::string _filter_max = "";

      bool _unique = false;
      int _rng_seed = 0;

      // More like system options.
      std::string _outdir = "output";

   };

#endif
